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

    ./Moore/run gaudirun.py hlt2_noUT_trackeff_test.py
"""

from Moore import options, run_moore
from Moore.tcks import dump_hlt2_configuration
from RecoConf.reconstruction_objects import reconstruction as reconstruction
from RecoConf.global_tools import stateProvider_with_simplified_geom
from RecoConf.decoders import default_ft_decoding_version
from Hlt2Conf.lines.trackeff.DiMuonTrackEfficiency import all_lines
import re

from Configurables import HiveDataBrokerSvc
HiveDataBrokerSvc().OutputLevel = 5

input_files = [
    #Channel Bs2JpsiPhi
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00146970/0000/00146970_00000003_1.xdigi',
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00146970/0000/00146970_00000006_1.xdigi',
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00146970/0000/00146970_00000007_1.xdigi',
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00146970/0000/00146970_00000008_1.xdigi',
    'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDIGI/00146970/0000/00146970_00000009_1.xdigi',
]
default_ft_decoding_version.global_bind(value=6)
options.dddb_tag = 'dddb-20210617'
options.conddb_tag = 'sim-20210617-vc-md100'
options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 0.3
options.evt_max = 100

options.output_file = 'hlt2_trackeff_noUT.dst'
options.output_type = 'ROOT'


def remove_lines(lines_dict, pattern_to_remove):
    filtered = {
        name: line
        for name, line in lines_dict.items()
        if re.search(pattern_to_remove, name) is None
    }
    return filtered


# Remove lines which need UT
muonut_to_remove = "Hlt2TrackEff_DiMuon_MuonUT.*"
downstream_to_remove = "Hlt2TrackEff_DiMuon_Downstream.*"

hlt2_lines = remove_lines(all_lines, muonut_to_remove)
trackeff_lines = remove_lines(hlt2_lines, downstream_to_remove)

print("Removed lines: ", all_lines.keys() - trackeff_lines.keys())

print("Number of HLT2 lines {}".format(len(trackeff_lines)))


def make_lines():
    return [builder() for builder in trackeff_lines.values()]


public_tools = [stateProvider_with_simplified_geom()]

with reconstruction.bind(from_file=False):
    config = run_moore(options, make_lines, public_tools)
dump_hlt2_configuration(config, "hlt2_trackeff_noUT.tck.json")
