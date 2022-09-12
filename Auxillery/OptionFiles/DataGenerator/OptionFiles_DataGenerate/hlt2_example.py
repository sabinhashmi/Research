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
    # MinBias 30000000
    # sim+std://MC/Upgrade/Beam7000GeV-Upgrade-MagDown-Nu7.6-25ns-Pythia8/Sim09c-Up02/Reco-Up01/30000000/LDST
    # 'root://xrootd.echo.stfc.ac.uk/lhcb:prod/lhcb/MC/Upgrade/LDST/00069155/0000/00069155_00000878_2.ldst'
    # D*-tagged D0 to KK, 27163002
    # sim+std://MC/Upgrade/Beam7000GeV-Upgrade-MagDown-Nu7.6-25ns-Pythia8/Sim09c-Up02/Reco-Up01/27163002/LDST
    #'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/LDST/00070317/0000/00070317_00000033_2.ldst'

    #'/data2/tracking/trackingGG/probki_danych/MC_Upgrade_Beam7000GeVUpgradeMagUpNu7.625nsPythia8_Sim10Up03OldP8Tuning_RecoUp02_11104142/00104527_00000032_2.xdst'
    '/data2/tracking/trackingGG/moore_build/MooreDev_v51r1/00104527_00000012_2.xdst'
#    '/data2/tracking/trackingGG/probki_danych/MC_Upgrade_Beam7000GeVUpgradeMagUpNu7.625nsPythia8_Sim10Up03OldP8Tuning_RecoUp02_11104142/00104527_00000012_2.xdst'


]

options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 4.2
# When running from Upgrade MC, must use the post-juggling locations of the raw
# event

options.evt_max = 1
options.simulation = True
options.data_type = 'Upgrade'
options.dddb_tag = 'dddb-20190726'
options.conddb_tag = 'sim-20190726-vc-mu100'
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
