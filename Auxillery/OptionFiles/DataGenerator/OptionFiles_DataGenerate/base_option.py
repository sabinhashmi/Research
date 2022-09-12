###############################################################################
# (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
"""Test persistreco when using real time reco. Produces hlt2_testpersistreco_realtime.dst

Run like any other options file:

    ./Moore/run gaudirun.py hlt2_persistreco_realtime.py
"""
import json

from Moore import options, run_moore
from Moore.config import HltLine
from Hlt2Conf.lines.charm.d0_to_hh import make_dzeros, make_charm_kaons
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.protoparticles import make_charged_protoparticles as make_charged_protoparticles_from
from RecoConf.reconstruction_objects import reconstruction, make_pvs, upfront_reconstruction
from RecoConf.hlt2_tracking import make_hlt2_tracks
from PyConf.Algorithms import PrGhostDumper, PrDownstreamDumper 
from PyConf.application import make_odin
from RecoConf.mc_checking import make_links_lhcbids_mcparticles_tracking_system, make_links_tracks_mcparticles, check_tracking_efficiency
from RecoConf.mc_checking_categories import get_mc_categories, get_hit_type_mask
from RecoConf.data_from_file import mc_unpackers
from RecoConf.hlt1_tracking import make_PrStoreUTHit_hits

################# Algorithms definitions #############################################
def ghost_dumper(hlt2_tracks):
    return PrGhostDumper(
            SeedTrackLocation = hlt2_tracks["Seed"]["v1"],
	    ODINLocation = make_odin())

def downstream_dumper(hlt2_tracks, seed_type = "Seed"):
    #with PrDownstreamDumper.bind(TrackLocation=hlt2_tracks["Seed"]["v1"]):
	    return PrDownstreamDumper(
        MCParticlesLocation = mc_unpackers()["MCParticles"],
        UTHitsLocation = make_PrStoreUTHit_hits(),
        ODINLocation = make_odin(),
#	LinkerLocation = make_links_lhcbids_mcparticles_tracking_system(),
	LinkerLocation = make_links_tracks_mcparticles(InputTracks=hlt2_tracks[seed_type]["v1"], LinksToLHCbIDs=make_links_lhcbids_mcparticles_tracking_system()),
	Tracks = hlt2_tracks[seed_type]["v1"])

def mc_checking(all_tracks, track_type="Downstream"):
	tracks = all_tracks[track_type]
	links_to_lhcbids = make_links_lhcbids_mcparticles_tracking_system()
	links_to_tracks = make_links_tracks_mcparticles( InputTracks=tracks, LinksToLHCbIDs=links_to_lhcbids)

	return check_tracking_efficiency(
			TrackType=track_type,
			InputTracks=tracks,
			LinksToTracks=links_to_tracks,
			LinksToLHCbIDs=links_to_lhcbids,
			MCCategories=get_mc_categories(track_type),
			HitTypesToCheck=get_hit_type_mask(track_type),
			)

def test_persistreco_line(name='Hlt2_test_persistrecoLine',
                          prescale=1,
                          persistreco=True):
    tracks = make_hlt2_tracks(light_reco=False)
    kaons = make_charm_kaons(pid_cut='PIDK > 15')
    seed_type = "Seed" # "Seed" / "Downstream"
    dzeros = make_dzeros(
        particles=kaons, descriptors=['D0 -> K- K+'], pvs=make_pvs())
    return HltLine(
        name=name,
        algs=upfront_reconstruction()\
			+[mc_checking(tracks)]\
			+[ghost_dumper(tracks)]\
			+[downstream_dumper(tracks, seed_type)]\
			,# + [dzeros],
        prescale=prescale,
        persistreco=persistreco)


###############################  Reco Options  #################################################


options.evt_max = 5

options.output_file = 'frozen_out.dst'
options.output_type = 'ROOT'
persist = False
remove_long_tracks = True


def make_lines():
    return [test_persistreco_line(persistreco=persist)]


from RecoConf.hlt2_tracking import make_hlt2_tracks
public_tools = [stateProvider_with_simplified_geom()]


with reconstruction.bind(from_file=False), \
	make_hlt2_tracks.bind(fast_reco=remove_long_tracks):
    run_moore(options, make_lines, public_tools)

# # temporarily save the PackedObjectLocations to get it back when reading in the file
# from Configurables import HltANNSvc
# with open('hlt2_persistreco_realtime.annsvc.pol.json', 'w') as outfile:
#     json.dump(
#         HltANNSvc().PackedObjectLocations, outfile, indent=4, sort_keys=True)

