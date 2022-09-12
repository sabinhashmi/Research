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
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000001_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000002_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000003_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000004_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000005_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000006_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000007_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000008_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000009_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000012_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000013_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000014_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000015_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000016_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000017_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000018_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000019_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000020_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000021_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000022_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000023_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000024_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000025_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000026_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000027_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000028_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000029_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000030_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000031_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000032_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000033_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000034_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000035_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000036_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000037_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000038_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000039_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000040_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000041_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000042_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000043_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000044_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000045_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000046_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000047_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000048_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000049_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000050_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000051_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000052_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000053_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000054_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000055_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000056_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000057_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000058_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000059_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000060_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000061_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000062_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000063_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000064_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000065_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000066_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000067_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000068_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000069_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000070_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000071_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000072_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000073_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000074_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000075_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000076_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000077_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000078_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000079_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000080_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000081_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000082_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000083_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000084_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000086_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000087_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000088_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000089_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000090_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000091_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/lambda/00102230_00000092_2.xdst',
]

options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 4.2
# When running from Upgrade MC, must use the post-juggling locations of the raw
# event

# options.evt_max = 10000
options.evt_max = 1
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
