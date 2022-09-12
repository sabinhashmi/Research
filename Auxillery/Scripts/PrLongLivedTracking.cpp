/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include <algorithm>
#include <bit>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "PrKernel/PrUTHitHandler.h"
#include "PrKernel/UTHitInfo.h"

#include "Event/PrDownstreamTracks.h"
#include "Event/PrSeedTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Kernel/ILHCbMagnetSvc.h"

#include "GaudiAlg/FixTESPath.h"
#include "GaudiAlg/Transformer.h"

#include "DetDesc/ConditionAccessorHolder.h" //FIXME

// from Gaudi
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <Math/CholeskyDecomp.h>
using ROOT::Math::CholeskyDecomp;

#include "UTDAQ/UTDAQHelper.h"
#include "UTDAQ/UTInfo.h"

// local
#include "PrDownTrack.h"

#include "weights/TMVA_PrLongLivedTracking_MLP.class.C"

#include "/home/hashmi/Files/TriggerPipeline/Auxillery/Model/GhostClassifier.cpp"

using DownTag = LHCb::Pr::Downstream::Tag;
using SeedTag = LHCb::Pr::Seeding::Tag;

/** @class PrLongLivedTracking PrLongLiveTracking.h
 *  Algorithm to reconstruct tracks with seed and UT, for daughters of long lived particles like Kshort, Lambda, ...
 *  Quick-and-dirty port from PatLongLivedTracking, more tuning needed
 *
 *  @author Michel De Cian and Adam Davis. Sascha Stahl and Olivier Callot for the original PrDownstream
 *  @date   2016-04-10
 *
 *  @2017-03-01: Christoph Hasse (adapt to future framework)
 *
 */
class PrLongLivedTracking
    : public Gaudi::Functional::Transformer<LHCb::Pr::Downstream::Tracks(
                                                const EventContext&, const LHCb::Pr::Seeding::Tracks&,
                                                const LHCb::Pr::UT::HitHandler&, const LHCb::UTDAQ::GeomCache& ),
                                            LHCb::DetDesc::usesConditions<LHCb::UTDAQ::GeomCache>> {
  using Track            = LHCb::Event::v2::Track;
  using SeedTracks       = LHCb::Pr::Seeding::Tracks;
  using DownstreamTracks = LHCb::Pr::Downstream::Tracks;

public:
  // - InputLocation: Input location of seed tracks
  // - OutputLocation: Output location of downstream tracks.
  // - ForwardLocation: Location of forward tracks, for hit rejection
  // - MatchLocation: Location of match tracks, for seed rejection
  PrLongLivedTracking( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"InputLocation", LHCb::TrackLocation::Seed}, KeyValue{"UTHits", UT::Info::HitLocation},
                      KeyValue{"GeometryInfo", "AlgorithmSpecific-" + name + "-UTGeometryInfo"}},
                     KeyValue{"OutputLocation", LHCb::TrackLocation::Downstream} ) {
    // flagging not supported yet
    // declareProperty( "ForwardLocation", m_ForwardTracks) ;
    // declareProperty( "MatchLocation", m_MatchTracks) ;
  }

  StatusCode initialize() override {
    return Transformer::initialize().andThen( [&] {
      // -- For MVA
      std::vector<std::string> mvaNames = {"chi2",      "seedChi2",    "p",     "pt",           "deltaP",
                                           "deviation", "initialChi2", "nHits", "highThresHits"};
      m_mvaReader                       = std::make_unique<NNForLLT>( mvaNames );

      addConditionDerivation<LHCb::UTDAQ::GeomCache( const DeUTDetector& )>( {DeUTDetLocation::UT},
                                                                             inputLocation<LHCb::UTDAQ::GeomCache>() );
    } );
  }

  DownstreamTracks operator()( const EventContext& evtCtx, const LHCb::Pr::Seeding::Tracks& InputTracks,
                               const LHCb::Pr::UT::HitHandler& hitHandler,
                               const LHCb::UTDAQ::GeomCache&   geometry ) const override {
    // create my state holding all needed mutable variables.
    std::array<Downstream::Hits, 4> preSelHits;
    Downstream::Hits                matchingXHits;
    Downstream::Hits                uHitsTemp;
    // -- track collections
    PrDownTracks goodXTracks;
    PrDownTracks goodXUTracks;
    PrDownTracks trackCandidates;

    matchingXHits.reserve( 64 );
    trackCandidates.reserve( 16 );
    goodXTracks.reserve( 8 );
    goodXUTracks.reserve( 8 );
    uHitsTemp.reserve( 64 );

    //==========================================================================
    // Get the output container
    //==========================================================================
    DownstreamTracks finalTracks{&InputTracks, Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
    finalTracks.reserve( InputTracks.size() );

    const double magScaleFactor = m_magFieldSvc->signedRelativeCurrent();

    bool magnetOff = std::abs( magScaleFactor ) > 1e-6 ? false : true;

    //==========================================================================
    // Main loop on tracks
    //==========================================================================
    for ( const auto& _tr : InputTracks.scalar() ) {
      const double model = catModel(_tr);
      if (model>0.5)continue;

      // -- simple Fisher discriminant to reject bad seed tracks
      // -- tune this!
      // const double fisher = evaluateFisher( tr );

      // if( fisher < m_seedCut ) continue;

      // -- Note: You want the state the furthest away from the magnet, as then
      // the straight-line approximation is the best
      const int          stateId = 2;
      const auto         state   = _tr.get<SeedTag::States>( stateId );
      Gaudi::TrackVector stateVector{state.x().cast(), state.y().cast(), state.tx().cast(), state.ty().cast(),
                                     _tr.qOverP( stateId ).cast()};
      PrDownTrack        refTrack( stateVector, state.z().cast(), m_zUT, m_zMagnetParams.value(), m_yParams.value(),
                            m_momentumParams.value(), magScaleFactor * ( -1 ) );

      if ( std::abs( refTrack.momentum() ) < 1400 ) continue;

      // -- Veto particles coming from the beam pipe.
      if ( insideBeampipe( refTrack ) ) continue;
      const double deltaP = refTrack.momentum() * refTrack.stateQoP() - 1.;

      // -- tune this!
      // -- check for compatible momentum
      if ( maxDeltaP( refTrack ) < fabs( deltaP ) ) {
        if ( !magnetOff ) continue;
      }

      // -- Get hits in UT around a first track estimate
      getPreSelection( refTrack, preSelHits, hitHandler, geometry );

      // -- Need at least 3 hits and at least 1 stereo and 1 x hit
      const int nXSelHits  = preSelHits[0].size() + preSelHits[3].size();
      const int nUVSelHits = preSelHits[1].size() + preSelHits[2].size();
      if ( 3 > nXSelHits + nUVSelHits ) continue;
      if ( 1 > nXSelHits ) continue;
      if ( 1 > nUVSelHits ) continue;

      trackCandidates.clear();

      //==============================================================
      // Try to find a candidate: X first, then UV.
      //==============================================================
      for ( int myPlane = 0; myPlane < 4; myPlane += 3 ) {
        const int otherPlane = ( myPlane == 0 ) ? 3 : 0;
        for ( auto& myHit : preSelHits[myPlane] ) {
          const double meanZ = myHit.z;
          const double posX  = myHit.x;

          PrDownTrack track( refTrack );

          // -- Create track estimate with one x hit
          const double slopeX = ( track.xMagnet() - posX ) / ( track.zMagnet() - meanZ );
          track.setSlopeX( slopeX );

          // -- Fit x projection
          findMatchingHits( track, preSelHits[otherPlane], matchingXHits );
          fitXProjection( track, myHit, matchingXHits, goodXTracks );

          // -- Loop over good x tracks
          for ( PrDownTrack& xTrack : goodXTracks ) {
            // -- Take all xTracks into account whose chi2 is close to the best
            if ( xTrack.chi2() - m_maxChi2DistXTracks >= goodXTracks[0].chi2() ) break;

            addUHits( xTrack, preSelHits[1], uHitsTemp, goodXUTracks );

            // -- Loop over good xu tracks
            for ( PrDownTrack& xuTrack : goodXUTracks ) {
              addVHits( xuTrack, preSelHits[2] );
              fitAndRemove<false>( xuTrack );
              if ( !acceptCandidate( xuTrack, magnetOff ) ) continue;
              trackCandidates.push_back( std::move( xuTrack ) );
            } // Loop over good xu tracks
          }   // Loop over good x tracks
        }
      }

      // Now we have all possible candidates, add overlap regions, fit again and
      // find best candidate

      PrDownTrack* bestCandidate    = nullptr;
      double       bestCandidateMVA = 0.05;

      for ( PrDownTrack& track : trackCandidates ) {

        addOverlapRegions( track, preSelHits );

        if ( track.chi2() > m_maxChi2 || insideBeampipe( track ) || track.hits().size() < 4 || track.hits().size() > 8 )
          continue;

        // -- Calculate the MVA
        double deviation = 0;

        for ( auto hit : track.hits() ) deviation += std::abs( track.distance( hit ) );

        // const int nHighThres = std::count_if( track.hits().begin(), track.hits().end(),
        //                                [](const Downstream::Hit& hit){ return hit.highThreshold(); });

        // -- This is a hack to test if we need this member in the UTHits.
        const int nHighThres = 3;

        const double initialChi2 = track.initialChi2();
        const double deltaP      = track.momentum() * track.stateQoP() - 1.;

        std::vector<double> vals = {track.chi2(),
                                    _tr.chi2PerDoF().cast(),
                                    std::abs( track.momentum() ),
                                    track.pt(),
                                    deltaP,
                                    deviation,
                                    initialChi2,
                                    static_cast<double>( track.hits().size() ),
                                    static_cast<double>( nHighThres )};

        const double mvaVal = m_mvaReader->GetMvaValue( vals );
        if ( bestCandidateMVA < mvaVal ) {
          bestCandidate    = &track;
          bestCandidateMVA = mvaVal;
        }
      }

      if ( bestCandidate ) { addTrack( *bestCandidate, _tr, finalTracks ); }
    } // Main loop on tracks

    m_downTrackCounter += finalTracks.size();

    return finalTracks;
  }

private:
  // Leave thes commented for now. These could possibly be used for as input for flagging, which is not supported
  // currently
  // DataObjectReadHandle<LHCb::Tracks> m_ForwardTracks  { this, "ForwardTracks" LHCb::TrackLocation::Forward  };
  // DataObjectReadHandle<LHCb::Tracks> m_MatchTracks    { this, "MatchTracks", LHCb::TrackLocation::Match   };
  //

  // - XPredTolConst: x-window for preselection is XPredTolConst/p + XPredTolOffset
  Gaudi::Property<double> m_xPredTolConst{this, "XPredTolConst", 200. * Gaudi::Units::mm* Gaudi::Units::GeV};
  // - XPredTolOffset: x-window for preselection is XPredTolConst/p + XPredTolOffset
  Gaudi::Property<double> m_xPredTolOffset{this, "XPredTolOffset", 6. * Gaudi::Units::mm};
  // - TolMatchConst: x-window for matching x hits is TolMatchConst/p + TolMatchOffset
  Gaudi::Property<double> m_tolMatchConst{this, "TolMatchConst", 20000.};
  //  - TolMatchOffset: x-window for matching x hits is TolMatchConst/p + TolMatchOffset
  Gaudi::Property<double> m_tolMatchOffset{this, "TolMatchOffset", 1.5 * Gaudi::Units::mm};
  // - TolUConst: window for u hits is TolUConst/p + TolUOffset
  Gaudi::Property<double> m_tolUConst{this, "TolUConst", 20000.0};
  // - TolUOffset: window for u hits is TolUConst/p + TolUOffset
  Gaudi::Property<double> m_tolUOffset{this, "TolUOffset", 2.5};
  // - TolVConst: window for v hits is TolVConst/p + TolVOffset
  Gaudi::Property<double> m_tolVConst{this, "TolVConst", 2000.0};
  // - TolVOffset: window for v hits is TolVConst/p + TolVOffset
  Gaudi::Property<double> m_tolVOffset{this, "TolVOffset", 0.5};
  // - MaxWindowSize: maximum window size for matching x hits
  Gaudi::Property<double> m_maxWindow{this, "MaxWindowSize", 10.0 * Gaudi::Units::mm};
  // - MaxChi2: Maximum chi2 for tracks with at least 4 hits
  Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 20.};
  // - MaxChi2ThreeHits: Maximum chi2 for tracks with 3 hits
  Gaudi::Property<double> m_maxChi2ThreeHits{this, "MaxChi2ThreeHits", 10.0};
  // - MinUTx: half-width  of beampipe rectangle
  Gaudi::Property<double> m_minUTx{this, "MinUTx", 25. * Gaudi::Units::mm};
  // - MinUTy: half-height of of beampipe rectangle
  Gaudi::Property<double> m_minUTy{this, "MinUTy", 25. * Gaudi::Units::mm};

  // Define parameters for MC09 field, zState = 9410
  // - ZMagnetParams: Parameters to determine the z-position of the magnet point. Tune with PrKsParams.
  Gaudi::Property<std::array<double, 7>> m_zMagnetParams{
      this, "ZMagnetParams", {5379.88, -2143.93, 366.124, 119074, -0.0100333, -0.146055, 1260.96}};
  // - MomentumParams: Parameters to determine the momentum. Tune with PrKsParams.
  Gaudi::Property<std::array<double, 3>> m_momentumParams{this, "MomentumParams", {1217.77, 454.598, 3353.39}};
  // - YParams: Parameters to determine the bending in y.  Tune with PrKsParams.
  Gaudi::Property<std::vector<double>> m_yParams{this, "YParams", {5., 2000.}};
  // - ZUT: z-position of middle of UT.
  Gaudi::Property<double> m_zUT{this, "ZUT", 2485. * Gaudi::Units::mm};
  // - ZUTa: z-position of first UT station
  Gaudi::Property<double> m_zUTa{this, "ZUTa", 2350. * Gaudi::Units::mm};
  // - MinPt: Minimum pT of the track
  Gaudi::Property<double> m_minPt{this, "MinPt", 0. * Gaudi::Units::MeV};
  // - MinMomentum: Minimum momentum of the track
  Gaudi::Property<double> m_minMomentum{this, "MinMomentum", 0. * Gaudi::Units::GeV};

  // -- Parameter to reject seed track which are likely ghosts
  // - FisherCut: Cut on Fisher-discriminant to reject bad seed tracks.
  Gaudi::Property<double> m_seedCut{this, "FisherCut", -1.0};

  // -- Parameters for the cut on deltaP (momentum estimate from Seeding and Downstream kink)
  // - MaxDeltaPConst: Window for deltaP is: MaxDeltaPConst/p + MaxDeltaPOffset
  Gaudi::Property<double> m_maxDeltaPConst{this, "MaxDeltaPConst", 0.0};
  // - MaxDeltaPOffset: Window for deltaP is: MaxDeltaPConst/p + MaxDeltaPOffset
  Gaudi::Property<double> m_maxDeltaPOffset{this, "MaxDeltaPOffset", 0.25};

  // -- Parameters for correcting the predicted position
  // - XCorrectionConst: Correction for x-position of search window in preselection is XCorrectionConst/p +
  //   XCorrestionOffset
  Gaudi::Property<double> m_xCorrectionConst{this, "XCorrectionConst", 23605.0};
  // - XCorrectionOffset: Correction for x-position of search window in preselection is XCorrectionConst/p +
  //   XCorrestionOffset
  Gaudi::Property<double> m_xCorrectionOffset{this, "XCorrestionOffset", 0.4};
  // - MaxXTracks: Maximum number of x-tracklets to process further
  Gaudi::Property<unsigned int> m_maxXTracks{this, "MaxXTracks", 2};
  // - MaxChi2DistXTracks: Maximum chi2 difference to x-tracklet with best chi2
  Gaudi::Property<double> m_maxChi2DistXTracks{this, "MaxChi2DistXTracks", 0.2};
  // - MaxXUTracks:  Maximum number of xu-tracklets to process further
  Gaudi::Property<unsigned int> m_maxXUTracks{this, "MaxXUTracks", 2};
  Gaudi::Property<double>       m_fitXProjChi2Offset{this, "FitXProjChi2Offset", 4.5};
  Gaudi::Property<double>       m_fitXProjChi2Const{this, "FitXProjChi2Const", 35000.0};

  // -- Tolerance for adding overlap hits
  // - OverlapTol: Tolerance for adding overlap hits
  Gaudi::Property<double> m_overlapTol{this, "OverlapTol", 2.0 * Gaudi::Units::mm};
  // - YTol: YTolerance for adding / removing hits.
  Gaudi::Property<double> m_yTol{this, "YTol", 2.0 * Gaudi::Units::mm};
  // Change this in order to remove hits and T-tracks used for longtracks.
  // RemoveAll configures that everything is removed.
  // If false only hits and T-tracks from good longtracks are removed.
  // The criterion for this is the Chi2 of the longtracks from the fit.
  // - RemoveUsed: Remove seed tracks and used UT hits (with chi2-cut on long track)?
  Gaudi::Property<bool> m_removeUsed{this, "RemoveUsed", false};
  // - RemoveAll: Remove seed tracks and used UT hits (withoug chi2-cut on long track)?
  Gaudi::Property<bool> m_removeAll{this, "RemoveAll", false};
  //  - LongChi2: Chi2-cut for the removal
  Gaudi::Property<double> m_longChi2{this, "LongChi2", 1.5};

  // void ttCoordCleanup() const;  ///< Tag already used coordinates

  // void prepareSeeds(const Tracks& inTracks, std::vector<Track*>& myInTracks)const; ///< Tag already used T-Seeds

  //=========================================================================
  //  Get the PreSelection of hits around a first track estimate
  //=========================================================================
  void getPreSelection( const PrDownTrack& track, std::array<Downstream::Hits, 4>& preSelHits,
                        const LHCb::Pr::UT::HitHandler& hitHandler, const LHCb::UTDAQ::GeomCache& geom ) const {
    // - Max Pt around 100 MeV for strange particle decay -> maximum displacement is in 1/p.
    double xPredTol = m_xPredTolOffset;

    // P dependance + overal tol.
    auto p = std::abs( track.momentum() );
    if ( p > 1e-6 ) xPredTol = m_xPredTolConst / p + m_xPredTolOffset;
    const double yTol = xPredTol / 2.0 + 7.5; // this is a little vague and not the final word

    // -- a correction turns out to be beneficial
    // -- maybe to compensate tracks not coming from 0/0 (?)
    const double correction = xPosCorrection( track );

    for ( auto& i : preSelHits ) i.clear();

    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();

    boost::container::small_vector<int, 9> sectors;

    for ( int iStation = 0; iStation < 2; ++iStation ) {
      if ( iStation == 1 && preSelHits[0].empty() && preSelHits[1].empty() ) return;
      for ( int iLayer = 0; iLayer < 2; ++iLayer ) {
        const int layer = 2 * iStation + iLayer;

        sectors.clear();

        const double zLayer = geom.layers[layer].z;
        const double yLayer = track.yAtZ( zLayer );
        const double xLayer = track.xAtZ( zLayer );
        LHCb::UTDAQ::findSectorsFullID( layer, xLayer, yLayer, xPredTol, yTol + 20.0, geom.layers[layer], sectors );
        std::sort( sectors.begin(), sectors.end() );

        int pp{-1};
        for ( auto& p : sectors ) {
          // sectors can be duplicated in the list, but they are ordered
          if ( p == pp ) continue;
          pp = p;

          const auto& sector  = hitHandler.indices( p );
          const auto& allHits = hitHandler.hits();

          if ( sector.first == sector.second ) continue;

          const auto myHs     = allHits.scalar();
          const auto firstHit = myHs[sector.first];

          const double dxDy     = firstHit.get<LHCb::Pr::UT::UTHitsTag::dxDy>().cast();
          const double zLayer   = firstHit.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>().cast();
          const double yPredLay = track.yAtZ( zLayer );
          const double xPredLay = track.xAtZ( zLayer );

          const double pos = xPredLay - correction;
          const auto   y   = yTrack + tyTr * zLayer;

          // this should sort of take the stereo angle and some tolerance into account.
          const double lowerBoundX = xPredLay - xPredTol - dxDy * yPredLay - 2.0;

          for ( auto itHit = sector.first; itHit < sector.second; ++itHit ) {

            const auto mH = myHs[itHit];

            const auto xAtYEq0 = mH.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast();

            const auto yBegin = mH.get<LHCb::Pr::UT::UTHitsTag::yBegin>().cast();
            const auto yEnd   = mH.get<LHCb::Pr::UT::UTHitsTag::yEnd>().cast();
            const auto yMin   = std::min( yBegin, yEnd );
            const auto yMax   = std::max( yBegin, yEnd );
            if ( xAtYEq0 < lowerBoundX || !( yMin - yTol <= yPredLay && yPredLay <= yMax + yTol ) ) continue;

            auto xx = mH.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast() + y * dxDy;
            if ( xPredTol < pos - xx ) continue; // go from -x to +x
            if ( xPredTol < xx - pos ) break;    // can break if we go out of the right bound

            preSelHits[layer].emplace_back( &allHits, itHit, xx, zLayer, fabs( xx - pos ) );
          }
        }
      }
    }

    std::sort( preSelHits[1].begin(), preSelHits[1].end(), Downstream::IncreaseByProj );
    std::sort( preSelHits[2].begin(), preSelHits[2].end(), Downstream::IncreaseByProj );
  }

  //=========================================================================
  //  Fit and remove the worst hit, as long as over tolerance
  //  Perform a chi2 fit to the track and remove outliers
  //=========================================================================
  template <bool onlyFit>
  void fitAndRemove( PrDownTrack& track ) const {
    if ( track.hits().size() < 2 ) {
      return; // no fit if single point only !
    }
    bool again = false;
    do {
      again = false;

      //== Fit, using the magnet point as constraint.
      double mat[6], rhs[3];
      mat[0] = track.weightXMag();
      mat[1] = 0.;
      mat[2] = 0.;
      mat[3] = 0.;
      mat[4] = 0.;
      mat[5] = 0.;
      rhs[0] = track.dxMagnet() * track.weightXMag();
      rhs[1] = 0.;
      rhs[2] = 0.;

      std::array<int, 4> differentPlanes = {0, 0, 0, 0};

      double yTrack = track.yAtZ( 0. );
      for ( auto& hit : track.hits() ) {
        if ( !onlyFit ) updateUTHitForTrackFast( hit, yTrack, track.slopeY() );

        const double dz   = 0.001 * ( hit.z - track.zMagnet() );
        const double dist = track.distance( hit );
        const double w    = hit.weight();
        const double t    = hit.sin();

        mat[0] += w;
        mat[1] += w * dz;
        mat[2] += w * dz * dz;
        mat[3] += w * t;
        mat[4] += w * dz * t;
        mat[5] += w * t * t;
        rhs[0] += w * dist;
        rhs[1] += w * dist * dz;
        rhs[2] += w * dist * t;

        // -- check how many different layers have fired
        differentPlanes[hit.planeCode()]++;
      }

      const unsigned int nDoF =
          std::count_if( differentPlanes.begin(), differentPlanes.end(), []( const int a ) { return a > 0; } );
      const int nbUV = differentPlanes[1] + differentPlanes[2];

      // -- solve the equation and update the parameters of the track
      CholeskyDecomp<double, 3> decomp( mat );
      if ( !decomp ) {
        track.setChi2( 1e42 );
        return;
      } else {
        decomp.Solve( rhs );
      }

      const double dx  = rhs[0];
      const double dsl = 0.001 * rhs[1];
      const double dy  = rhs[2];

      if ( 4 > nbUV ) track.updateX( dx, dsl );
      track.updateY( dy );

      //== Remove worst hit and retry, if too far.
      double chi2 = track.initialChi2();

      double                      maxDist = -1.;
      PrDownTrack::Hits::iterator worst;

      yTrack = track.yAtZ( 0. );
      for ( auto itH = track.hits().begin(); itH != track.hits().end(); ++itH ) {
        Downstream::Hit& hit = *itH;
        updateUTHitForTrackFast( hit, yTrack, track.slopeY() );

        if ( !onlyFit ) {
          const double yTrackAtZ = track.yAtZ( hit.z );
          if ( !hit.isYCompatible( yTrackAtZ, m_yTol ) ) {
            track.hits().erase( itH );
            if ( 2 < track.hits().size() ) again = true;
            break;
          }
        }

        const double dist = std::abs( track.distance( hit ) );
        hit.projection    = dist;
        chi2 += dist * dist * hit.weight();
        // -- only flag this hit as removable if it is not alone in a plane or there are 4 planes that fired
        if ( !onlyFit && maxDist < dist && ( 1 < differentPlanes[hit.planeCode()] || nDoF == track.hits().size() ) ) {
          maxDist = dist;
          worst   = itH;
        }
      }

      if ( again ) continue;

      if ( 2 < track.hits().size() ) chi2 /= ( track.hits().size() - 2 );
      track.setChi2( chi2 );

      if ( onlyFit ) { return; }

      if ( m_maxChi2 < chi2 && track.hits().size() > 3 && maxDist > 0 ) {
        track.hits().erase( worst );
        again = true;
      }
    } while ( again );
  }

  //=========================================================================
  //  Collect the hits in the other x layer
  //=========================================================================
  void findMatchingHits( const PrDownTrack& track, const Downstream::Hits& preSelHits,
                         Downstream::Hits& matchingXHits ) const {
    matchingXHits.clear();
    if ( preSelHits.empty() ) return;

    const double tol =
        std::min( m_maxWindow.value(), m_tolMatchOffset + m_tolMatchConst * std::abs( track.stateQoP() ) );
    const double xPred = track.xAtZ( preSelHits.front().z );

    for ( auto& hit : preSelHits ) {
      const double adist = std::abs( hit.x - xPred );
      if ( adist <= tol ) { matchingXHits.push_back( hit ); }
    }

    std::sort( matchingXHits.begin(), matchingXHits.end(),
               [xPred]( const Downstream::Hit& lhs, const Downstream::Hit& rhs ) {
                 return std::abs( lhs.x - xPred ) < std::abs( rhs.x - xPred );
               } );
  }

  //=========================================================================
  //  Store Track
  //=========================================================================
  template <typename ProxyType>
  void addTrack( const PrDownTrack& track, const ProxyType& seed, LHCb::Pr::Downstream::Tracks& tracks ) const {
    auto newTrack = tracks.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

    newTrack.field<DownTag::trackSeed>().set( seed.offset() );

    // Store UT hits
    int i = 0;
    for ( const auto& hit : track.hits() ) {
      LHCb::Event::lhcbid_v<SIMDWrapper::scalar::types> id( static_cast<unsigned int>( hit.lhcbID() ) );
      newTrack.field<DownTag::lhcbIDs>( i ).setLHCbID( id );
      newTrack.field<DownTag::ut_indices>( i ).set( hit.hit );
      i++;
    }
    newTrack.field<DownTag::nUTHits>().set( i );

    // Copy seed hits
    for ( i = 0; i < seed.nHits().cast(); i++ ) {
      auto id = seed.template field<LHCb::Pr::Seeding::Tag::lhcbIDs>( i ).lhcbID();
      // newTrack.field<DownTag::lhcbIDs>( track.hits().size() + i ).setLHCbID( seed.lhcbID( i ) );
      newTrack.field<DownTag::lhcbIDs>( track.hits().size() + i ).setLHCbID( id );
      newTrack.field<DownTag::ft_indices>( i ).set( seed.ft_index( i ) );
    }
    newTrack.field<DownTag::nFTHits>().set( seed.nHits() );

    // Create a state at zUTa
    newTrack.field<DownTag::State>().setPosition( static_cast<float>( track.xAtZ( m_zUTa ) ),
                                                  static_cast<float>( track.yAtZ( m_zUTa ) ),
                                                  static_cast<float>( m_zUTa ) );
    newTrack.field<DownTag::State>().setDirection( static_cast<float>( track.slopeX() ),
                                                   static_cast<float>( track.slopeY() ) );
    newTrack.field<DownTag::State>().setQOverP( static_cast<float>( 1.0 / track.momentum() ) );
  }

  //=========================================================================
  //  Add the U hits.
  //=========================================================================
  void addUHits( const PrDownTrack& track, Downstream::Hits& preSelHits, Downstream::Hits& uHitsTemp,
                 PrDownTracks& goodXUTracks ) const {
    goodXUTracks.clear();

    if ( preSelHits.empty() ) { return; }

    uHitsTemp.clear();

    const double tol = m_tolUOffset + m_tolUConst / std::abs( track.momentum() );

    // -- these numbers are a little arbitrary
    double minChi2 = ( track.hits().size() == 1 ) ? 800 : 300;

    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();

    // -- first select all hits, and then
    // -- accept until over a tolerance
    for ( auto& hit : preSelHits ) {
      updateUTHitForTrackFast( hit, yTrack, tyTr );

      const double dist = std::abs( track.distance( hit ) );
      if ( dist > tol ) continue;
      hit.projection = dist;
      uHitsTemp.push_back( hit );
    }

    std::sort( uHitsTemp.begin(), uHitsTemp.end(), Downstream::IncreaseByProj );

    for ( auto& hit : uHitsTemp ) {
      PrDownTrack greatTrack( track );

      greatTrack.hits().push_back( hit );
      fitAndRemove<true>( greatTrack );

      // -- it's sorted
      if ( greatTrack.chi2() > minChi2 ) break;
      goodXUTracks.push_back( std::move( greatTrack ) );
      if ( goodXUTracks.size() >= m_maxXUTracks ) { break; }
    }
  }

  //=========================================================================
  //  Add the V hits. Take the one which has the best chi2
  //=========================================================================
  void addVHits( PrDownTrack& track, Downstream::Hits& preSelHits ) const {
    if ( preSelHits.empty() ) { return; }

    auto   p   = std::abs( track.momentum() );
    double tol = ( track.hits().size() == 2 ) ? m_tolUOffset + m_tolUConst / p : m_tolVOffset + m_tolVConst / p;

    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();

    double minChi2 = 10000;

    Downstream::Hit* bestHit = nullptr;
    for ( auto& hit : preSelHits ) {
      updateUTHitForTrackFast( hit, yTrack, tyTr );
      const double adist = std::abs( track.distance( hit ) );

      if ( adist < tol ) {
        hit.projection = adist;
        track.hits().push_back( hit );
        fitAndRemove<true>( track );
        track.hits().pop_back();

        if ( track.chi2() < minChi2 ) {
          bestHit = &hit;
          minChi2 = track.chi2();
        }
      }
    }

    if ( bestHit ) track.hits().push_back( *bestHit );

    track.sortFinalHits();
  }

  // void tagUsedUT( const Track* tr ) const; ///< Tag hits that were already used elsewhere

  //=============================================================================
  // Fit the projection in the zx plane, one hit in each x layer
  //=============================================================================
  void fitXProjection( const PrDownTrack& track, Downstream::Hit& firstHit, Downstream::Hits& matchingXHits,
                       PrDownTracks& goodXTracks ) const {
    goodXTracks.clear();

    const double maxChi2 = m_fitXProjChi2Offset + m_fitXProjChi2Const / std::abs( track.momentum() );

    // Catch if there is no second hit in other station
    for ( const auto& hit : matchingXHits ) {
      PrDownTrack tr( track );
      xFit( tr, firstHit, hit );

      if ( tr.chi2() > maxChi2 ) { break; }
      tr.hits().push_back( firstHit );
      tr.hits().push_back( hit );

      goodXTracks.push_back( std::move( tr ) );
      if ( goodXTracks.size() >= m_maxXTracks ) { break; }
    }

    std::sort( goodXTracks.begin(), goodXTracks.end(),
               []( const PrDownTrack& a, const PrDownTrack& b ) { return a.chi2() < b.chi2(); } );
  }

  //=========================================================================
  //  Check if the new candidate is better than the old one
  //=========================================================================
  bool acceptCandidate( const PrDownTrack& track, bool magnetOff ) const {
    const int nbMeasureOK = track.hits().size();

    //== Enough measures to have Chi2/ndof.
    if ( nbMeasureOK < 3 ) { return false; }

    // -- use a tighter chi2 for 3 hit tracks
    // -- as they are more likely to be ghosts
    const double maxChi2 = ( nbMeasureOK == 3 ) ? m_maxChi2ThreeHits : m_maxChi2;

    //== Good enough Chi2/ndof
    if ( maxChi2 < track.chi2() ) { return false; }

    //== Compatible momentum
    const double p      = track.momentum();
    const double deltaP = p * track.stateQoP() - 1.;

    if ( maxDeltaP( track ) < fabs( deltaP ) && !magnetOff ) { return false; }
    if ( std::abs( p ) < m_minMomentum || track.pt() < m_minPt ) { return false; }

    return true;
  }

  //=============================================================================
  // This is needed for tracks which have more than one x hit in one layer
  // Maybe we could make this smarter and do it for every track and add the 'second best'
  // this, such that we do not need to loop over them again
  //=============================================================================
  void addOverlapRegions( PrDownTrack& track, std::array<Downstream::Hits, 4>& preSelHits ) const {
    bool hitAdded = false;

    const double yTrack      = track.yAtZ( 0. );
    const int    trackNbHits = track.hits().size();

    for ( int i = 0; i < 4; ++i ) {
      for ( auto& hit : preSelHits[i] ) {

        updateUTHitForTrackFast( hit, yTrack, track.slopeY() );

        if ( m_overlapTol > std::abs( track.distance( hit ) ) ) {
          double yTrack = track.yAtZ( hit.z );

          if ( !hit.isYCompatible( yTrack, m_yTol ) ) continue;

          // -- check that z-position is different
          bool addHit = true;
          int  j      = 0;
          for ( const auto& trackHit : track.hits() ) {
            if ( j >= trackNbHits ) break;
            j++;
            if ( trackHit.planeCode() != hit.planeCode() ) continue;
            // -- the displacement in z between overlap modules is larger than 1mm
            if ( std::abs( hit.z - trackHit.z ) < 1.0 ) {
              addHit = false;
              break;
            }
          }

          if ( addHit ) {
            track.hits().push_back( hit );
            hitAdded = true;
          }
        }
      }
    }

    if ( hitAdded ) {
      track.sortFinalHits();
      fitAndRemove<false>( track );
    }
  }

  //=============================================================================
  // Evaluate the Fisher discriminant for a preselection of seed tracks
  //=============================================================================
  template <typename ProxyType_> 
  double catModel( const ProxyType_& seed ) const{
    

    int          stateId = 2;
    LHCb::State state = seed.getLHCbState(stateId);

    

    double chi2PerDoF = seed.chi2PerDoF().cast();

    double nlhcbID = seed.nHits().cast();

    double p = state.p();

    double pt = state.pt();
    
    double phi = state.momentum().phi();

    double position_x = state.x();
    
    double position_y = state.y();
    
    double position_r = sqrt(pow(position_x,2)+pow(position_y,2));

    double tx = state.tx();
    double ty = state.ty();
    
    
    double pseudoRapidity = state.momentum().eta();
    
    // double qOverP = seed.qOverP().cast();
    

    std::vector<double> parameters = {chi2PerDoF, nlhcbID, p, phi , position_r, position_x, position_y, pt, tx, ty, pseudoRapidity};




    auto model_out = ApplyCatboostModel(parameters);
    auto probability = exp(model_out) / (1+exp(model_out));

    // std::fstream dump_file;
    // dump_file.open("/home/hashmi/Files/StandaloneAnalysis/ThresholdCheck_CompleteData/Dumps/Output.txt",std::fstream::in | std::fstream::out | std::fstream::app);

    // dump_file << p << "," << chi2PerDoF << "," << nlhcbID << ","  << phi << ","  << position_r << "," << position_x << "," << position_y << "," << pt << "," << tx << ","<< ty << "," << pseudoRapidity << "," << qOverP << "," << probability <<std::endl;

    // dump_file << model_out << "," << probability <<std::endl;

    // dump_file.close();

    return probability;
  }

  //=========================================================================
  //  Fit hits in x layers, using the magnet point as constraint.
  //=========================================================================
  void xFit( PrDownTrack& track, const Downstream::Hit& hit1, const Downstream::Hit& hit2 ) const {
    double mat[3], rhs[2];

    const auto w1 = hit1.weight();
    const auto w2 = hit2.weight();

    auto d1 = track.distance( hit1 );
    auto d2 = track.distance( hit2 );

    const auto dz1 = hit1.z - track.zMagnet();
    const auto dz2 = hit2.z - track.zMagnet();

    mat[0] = track.weightXMag() + w1 + w2;
    mat[1] = w1 * dz1 + w2 * dz2;
    mat[2] = w1 * dz1 * dz1 + w2 * dz2 * dz2;

    rhs[0] = w1 * d1 + w2 * d2;
    rhs[1] = w1 * d1 * dz1 + w2 * d2 * dz2;

    // Solve linear system
    const double det = mat[0] * mat[2] - mat[1] * mat[1];
    const double dx  = ( mat[2] * rhs[0] - mat[1] * rhs[1] ) / det;
    const double dsl = ( mat[0] * rhs[1] - mat[1] * rhs[0] ) / det;

    track.updateX( dx, dsl );

    d1 = track.distance( hit1 );
    d2 = track.distance( hit2 );

    const double chi2 = track.initialChi2() + w1 * d1 * d1 + w2 * d2 * d2;
    track.setChi2( chi2 );
  }

  // -- UT hits don't need the correction for the differnce in the reference frame, this should be a bit faster
  // -- than the original implementation
  void updateUTHitForTrackFast( Downstream::Hit& hit, const double y0, const double dyDz ) const {
    const auto y = ( y0 + dyDz * hit.zAtYEq0() );
    hit.x        = hit.xAt( y );
  }

  /// Does this track point inside the beampipe?
  bool insideBeampipe( const PrDownTrack& track ) const {
    return ( m_minUTx > fabs( track.xAtZ( m_zUTa ) ) ) && ( m_minUTy > fabs( track.yAtZ( m_zUTa ) ) );
  }

  /// Helper to evaluate the Fisher discriminant
  /* double getFisher(const std::array<double,5> vals) {
    double                c_fishConst        = -1.69860581797;
    std::array<double, 5> c_fishCoefficients = {-0.241020410138, 3.03197732663e-07, -1.14400162824e-05,
                                                 0.126857153245, 0.122359738469};
    double fishVal = c_fishConst;
    for(int i = 0; i < 5; i++) fishVal += vals[i] * c_fishCoefficients[i];
    return fishVal;
  } */

  /// Helper to evaluate the maximum discrepancy between momentum from kink and curvature in T-stations
  double maxDeltaP( const PrDownTrack& track ) const {
    return m_maxDeltaPConst / std::abs( track.momentum() ) + m_maxDeltaPOffset;
  }

  /// Helper to evaluate the correction to the x position in the UT
  double xPosCorrection( const PrDownTrack& track ) const {
    auto p = track.momentum();
    return std::copysign( m_xCorrectionOffset + m_xCorrectionConst / std::abs( p ), p );
  }

  // -- counters
  mutable Gaudi::Accumulators::SummingCounter<> m_downTrackCounter{this, "#Downstream tracks made"};
  mutable Gaudi::Accumulators::SummingCounter<> m_utHitsCounter{this, "#UT hits added"};

  // -- pointers to tools
  ServiceHandle<ILHCbMagnetSvc>            m_magFieldSvc{this, "MagneticField", "MagneticFieldSvc"};
  std::unique_ptr<IClassifierReaderForLLT> m_mvaReader;
};

DECLARE_COMPONENT( PrLongLivedTracking )
