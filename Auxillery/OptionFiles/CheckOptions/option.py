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
"""Options for running HLT2 lines.

Run like any other options file:

    ./Moore/run gaudirun.py hlt2_example.py
"""
from Moore import options, run_moore

#from Hlt1Conf.lines.high_pt_muon import all_lines
#from RecoConf.global_tools import stateProvider_with_simplified_geom
# TODO stateProvider_with_simplified_geom must go away from option files

input_files = [
# 'LFN:/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000005_2.xdst'
# 'sim+std://MC/Upgrade/Beam7000GeV-Upgrade-MagDown-Nu7.6-25ns-Pythia8/Sim10-Up03-OldP8Tuning/Reco-Up02/11104142/XDST'
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/LDST/00070317/0000/00070317_00000033_2.ldst'
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000005_2.xdst'
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000005_2.xdst'
]

options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 4.3
options.n_threads= 1
# When running from Upgrade MC, must use the post-juggling locations of the raw
# event

# options.evt_max = 10000
options.evt_max = 10
options.simulation = True
options.data_type = 'Upgrade'
options.dddb_tag = 'dddb-20190223'
options.conddb_tag = 'sim-20180530-vc-md100'
#options.print_freq = 1

from RecoConf.hlt1_tracking import default_ft_decoding_version
default_ft_decoding_version.global_bind(value=6)

#from Configurables import VPClus, VeloClusterTrackingSIMD, DigiConf, FTRawBankDecoder
#VPClus("VPClustering").RawEventLocation = "/Event/Velo/RawEvent"
#VeloClusterTrackingSIMD("VeloClusterTracking").RawEventLocation = "/Event/Velo/RawEvent"
#FTRawBankDecoder("createFtClusters").DecodingVersion = 6

#DigiConf().EnableUnpack = True


#def make_lines():
#    return [builder() for builder in all_lines.values()]


#public_tools = [stateProvider_with_simplified_geom()]
#run_reconstruction(options, make_lines, public_tools)
#run_moore(options, make_lines, public_tools)



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


from Moore import  run_reconstruction
from RecoConf.standalone import *
from RecoConf.hlt2_tracking import make_hlt2_tracks
from PyConf.Algorithms import PrLongLivedTracking
rec = standalone_hlt2_full_track_reco


with rec.bind(do_mc_checking=True):
	with make_hlt2_tracks.bind(fast_reco=True):
		run_reconstruction(options, rec)


