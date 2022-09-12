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

from Moore import options

input_files = [
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000021_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000043_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000026_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000029_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000042_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000025_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000049_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000054_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000058_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000063_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000055_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000076_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000086_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000088_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000007_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000031_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000027_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000036_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000047_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000066_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000065_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000094_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000005_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000019_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000015_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000035_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000039_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000060_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000052_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000062_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000071_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000099_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000004_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000008_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000003_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000016_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000030_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000034_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000024_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000046_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000053_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000057_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000056_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000061_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000084_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000083_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000085_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000095_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000091_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000006_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000009_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000001_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000040_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000028_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000044_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000048_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000064_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000069_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000080_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000082_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000090_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000010_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000020_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000023_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000011_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000017_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000033_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000041_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000070_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000077_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000078_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000081_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000002_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000013_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000018_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000014_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000012_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000032_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000089_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000022_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000037_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000038_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000045_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000051_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000072_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000059_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000073_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000067_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000068_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000074_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000075_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000079_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000087_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000096_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000097_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000093_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000092_2.xdst',
'/data2/tracking/trackingGG/frozen_build/data/ks/00104521_00000098_2.xdst',

]
options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 4.3

options.simulation = True
options.data_type = 'Upgrade'
options.dddb_tag = 'dddb-20190223'
options.conddb_tag = 'sim-20180530-vc-md100'


from RecoConf.hlt1_tracking import default_ft_decoding_version
default_ft_decoding_version.global_bind(value=6)

#-- GAUDI jobOptions generated on Wed Jan 26 19:23:35 2022
#-- Contains event types : 
#--   11104142 - 98 files - 103521 events - 396.46 GBytes

#--  Extra information about the data processing phases:

#--  Processing Pass: '/Sim10-Up03-OldP8Tuning/Reco-Up02' 

#--  StepId : 139214 
#--  StepName : Digi15-Upgrade for Upgrade studies with spillover - 2017 Baseline NoRichSpillover - xdigi 
#--  ApplicationName : Boole 
#--  ApplicationVersion : v40r4 
#--  OptionFiles : $APPCONFIGOPTS/Boole/Default.py;$APPCONFIGOPTS/Boole/Boole-Upgrade-Baseline-20150522.py;$APPCONFIGOPTS/Boole/EnableSpillover.py;$APPCONFIGOPTS/Boole/Upgrade-RichMaPMT-NoSpilloverDigi.py;$APPCONFIGOPTS/Boole/xdigi.py 
#--  DDDB : fromPreviousStep 
#--  CONDDB : fromPreviousStep 
#--  ExtraPackages : AppConfig.v3r376 
#--  Visible : N 

#--  Processing Pass: '/Sim10-Up03-OldP8Tuning/Reco-Up02' 

#--  StepId : 140462 
#--  StepName : Reco-Up02 for Upgrade studies - 2017 Baseline Geometry - XDST 
#--  ApplicationName : Brunel 
#--  ApplicationVersion : v60r6p1 
#--  OptionFiles : $APPCONFIGOPTS/Brunel/MC-WithTruth.py;$APPCONFIGOPTS/Brunel/Brunel-Upgrade-Baseline-20150522.py;$APPCONFIGOPTS/Brunel/patchUpgrade1.py;$APPCONFIGOPTS/Brunel/xdst.py 
#--  DDDB : fromPreviousStep 
#--  CONDDB : fromPreviousStep 
#--  ExtraPackages : AppConfig.v3r392 
#--  Visible : Y 