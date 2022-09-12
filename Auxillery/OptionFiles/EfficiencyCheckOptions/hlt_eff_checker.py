#!/usr/bin/env python
###############################################################################
# (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from __future__ import absolute_import, print_function
import argparse
import os
import yaml
from datetime import datetime
from HltEfficiencyChecker.wizard import main

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=
        "Script that extracts and plots efficiencies and rates from the high level trigger."
    )

    parser.add_argument("config", help="YAML configuration file")

    parser.add_argument(
        "-o",
        "--output",
        default="checker",
        help="Output directory prefix. See also --output-suffix.")

    parser.add_argument(
        "-s",
        "--output-suffix",
        default="-{now:%Y%m%d-%H%M%S}",
        help="Output directory suffix")

    parser.add_argument(
        "--force",
        action="store_true",
        help="Do not fail when output directory exists")

    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Only print the commands needed to run from stack/ directory.")

    parser.add_argument(
        "--ignore-broken-inputs",
        action="store_true",
        help=
        "Ignore Gaudi::ReturnCode::FailInput errors in check_call. Tupling already skips these files, so this option allows the rest of HltEfficiencyChecker to proceed after skipping a broken input file."
    )

    args = parser.parse_args()

    args.output_suffix = args.output_suffix.format(now=datetime.now())
    cwd = args.output + args.output_suffix if args.output else ""

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    run_command = os.getcwd() + '/MooreAnalysis/run'
    if cwd:
        if not os.path.isdir(cwd):
            os.mkdir(cwd)
        elif not args.force:
            parser.error(
                "Output directory {} exists. Use --force to use anyway".format(
                    cwd))
        if args.dry_run:
            print("The commands to run are... ")
        print("cd {!r}".format(cwd))
        os.chdir(cwd)

    main(
        config,
        run_command,
        dry_run=args.dry_run,
        ignore_broken_inputs=args.ignore_broken_inputs)
