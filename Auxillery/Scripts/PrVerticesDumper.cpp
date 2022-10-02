#include "Event/ODIN.h"
#include "Event/RecVertex_v2.h"

#include "GaudiAlg/Consumer.h"
#include "boost/filesystem.hpp"
#include <utility>

#include "TROOT.h"
#include <TFile.h>
#include <TString.h>
#include <TTree.h>


class PrVerticesDumper
	: public Gaudi::Functional::Consumer<void( const LHCb::Event::v2::RecVertices&, const LHCb::ODIN&) >{
		public:
			PrVerticesDumper( const std::string& name, ISvcLocator* pSvcLocator );
			StatusCode initialize() override;
			void operator()( const LHCb::Event::v2::RecVertices&, const LHCb::ODIN& )const override;

		private:
			Gaudi::Property<std::string> m_output_directory{this, "OutputDirectory", ".", "Directory to store output nTuples"};
	
	};
DECLARE_COMPONENT( PrVerticesDumper )

namespace {
	namespace fs = boost::filesystem;
}

PrVerticesDumper::PrVerticesDumper( const std::string& name, ISvcLocator* pSvcLocator )
	: Consumer( name, pSvcLocator,
			{KeyValue{"VerticesLocation", LHCb::Event::v2::RecVertexLocation::Primary},
			KeyValue{"ODINLocation", LHCb::ODINLocation::Default}} ) {}
StatusCode PrVerticesDumper::initialize() {
	StatusCode sc = Consumer::initialize();
	if (sc.isFailure() ) return sc;

	auto dir = fs::path{m_output_directory.value()};

	if ( !fs::exists(dir)){
		boost::system::error_code ec;
		bool success = fs::create_directories(dir,ec);
		success &= !ec;
		if(!success){
			error() << "Failed to create directory " << dir.string() << ": " << ec.message() << endmsg;
			return StatusCode::FAILURE;
		}
	}
	return sc;
}

void PrVerticesDumper::operator()( const LHCb::Event::v2::RecVertices& Vertices, const LHCb::ODIN& odin)const {
	std::string filename = ( "VerticesDumper_runNb_"+std::to_string(odin.runNumber()) + "_evtNb_" + std::to_string(odin.eventNumber())+".root" );
	TFile file( filename.c_str(), "UPDATE");

	int eventID = 0;
	double vertex_x = 0.0;
	double vertex_y = 0.0;
	double vertex_z = 0.0;

	auto tree = TTree("VertexData", "VertexData");
	tree.Branch( "vertex_x", &vertex_x);
	tree.Branch( "vertex_y", &vertex_y);
	tree.Branch( "vertex_z", &vertex_z);
	tree.Branch( "eventID", &eventID);

	eventID = odin.eventNumber();
	for (auto const& vertex: Vertices){
		vertex_x = vertex.position().x();
		vertex_y = vertex.position().y();
		vertex_z = vertex.position().z();
		tree.Fill();
	}

	file.Write();
	file.Close();


}


