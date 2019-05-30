#!/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import optparse
import shutil
import time
import seqnum

# Import Broadband modules
from install_cfg import InstallCfg
import bband_utils

# Enumerated types
class Methods(object):
    """
    Defines available models on the platform
    """
    gp, ucsb, sdsu, exsim, song, irikura1, irikura2, csm = range(8)
    labels = ["GP", "UCSB", "SDSU", "EXSIM", "SONG", "IRIKURA_RECIPE_M1", "IRIKURA_RECIPE_M2"]
    options = ["gp", "ucsb", "sdsu", "exsim", "song", "irikura1", "irikura2"]
    #labels = ["GP", "UCSB", "SDSU", "EXSIM", "SONG", "IRIKURA_RECIPE_M1", "IRIKURA_RECIPE_M2", "CSM"]
    #options = ["gp", "ucsb", "sdsu", "exsim", "song", "irikura1", "irikura2", "csm"]

class GenAcceptTests(object):
    def __init__(self, resume=True):
        install = InstallCfg.getInstance()
        self.resume = resume
        self.resume_list = []

        # Read checkpoint file
        if self.resume == True:
            resume_file = os.path.join(install.A_OUT_LOG_DIR,
                                       "gen_resume.txt")
            if os.path.exists(resume_file):
                resume_fp = open(resume_file, 'r')
                self.resume_list = resume_fp.read().splitlines()
                resume_fp.close()
            else:
                self.resume = False

        # Setup paths
        self.input_dir = os.path.join(install.A_TEST_REF_DIR,
                                      "accept_inputs")
        self.ref_dir = os.path.join(install.A_TEST_REF_DIR,
                                    "accept_refs")

    def generate_option_list(self):
        """
        This function creates the option files for all test cases
        """
        install = InstallCfg.getInstance()

        numvalid = 0
        numuser = 0
        optfiles = {}

        # Make sure run directory exists
        cmd = "mkdir -p %s" % (install.A_USER_DATA_DIR)
        bband_utils.runprog(cmd)

        # Clear out run directory
        user_temp = os.path.join(install.A_USER_DATA_DIR, "tmp")
        if os.path.exists(user_temp):
            shutil.rmtree(user_temp)
            filelist = os.listdir(install.A_USER_DATA_DIR)
            if not os.path.exists(user_temp):
                os.mkdir(user_temp)
            for file_entry in filelist:
                shutil.move(os.path.join(install.A_USER_DATA_DIR, file_entry),
                            os.path.join(user_temp, file_entry))

        basedir = os.path.join(install.A_TEST_REF_DIR, "accept_inputs")

        # Copy input files to run directory
        shutil.copy2(os.path.join(basedir, "northridge_3_sta.stl"),
                     "%s" % install.A_USER_DATA_DIR)
        shutil.copy2(os.path.join(basedir, "northridge_eq_gp.src"),
                     "%s" % install.A_USER_DATA_DIR)
        shutil.copy2(os.path.join(basedir, "northridge_eq_ucsb.src"),
                     "%s" % install.A_USER_DATA_DIR)
        shutil.copy2(os.path.join(basedir, "northridge_eq_song.src"),
                     "%s" % install.A_USER_DATA_DIR)

        # Ensure input dir and ref dir exist
        cmd = "mkdir -p %s" % (self.input_dir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (self.ref_dir)
        bband_utils.runprog(cmd)

        # Get sorted list of station files for menu selections
        files = glob.glob("%s%s*.stl" % (install.A_USER_DATA_DIR, os.sep))
        stafiles = []
        for file_entry in files:
            stafiles.append(os.path.basename(file_entry))
        stafiles.sort()

        # Validation simulations
        mode = "valid-northridge"
        for method in xrange(0, len(Methods.labels)):
            optfile = "%s-%s.txt" % (mode, Methods.labels[method])
            print("Generating %s" % (optfile))
            opts = []
            # Select validation run
            opts.append('y')
            # Select the validation event
            opts.append("NR")
            # Select method
            opts.append(Methods.options[method])
            # We don't need a custom source file
            opts.append('n')
            # For GP, UCSB, and SDSU, SONG, Irikura1, Irikura2
            # we want to run the rupture generator
            if method in [Methods.gp, Methods.ucsb, Methods.sdsu, Methods.song, Methods.irikura1, Methods.irikura2]:
                opts.append('y')
            opts.append('2')
            opts.append('1')
            opts.append("%d" %
                        (stafiles.index("northridge_3_sta.stl") + 1))
            if method == Methods.exsim:
                # No custom EXSIM template file
                opts.append('n')
            if method != Methods.csm:
                # Run site response (CSM does not ask this question)
                opts.append('y')
            # Skip plots
            opts.append('n')
            opts.append('n')
            # Skip GMPE Comparison
            opts.append('n')
            # Yes to GoF
            opts.append('y')
            # GP GoF
            opts.append('1')
            # No to additional metrics
            opts.append('n')
            print("\t %s" % (str(opts)))
            numvalid = numvalid + 1
            optfiles[optfile] = opts

        # User simulations
        mode = "user"
        for method in xrange(0, len(Methods.labels)):
            optfile = "%s-%s.txt" % (mode, Methods.labels[method])
            print("Generating %s" % (optfile))
            opts = []
            opts.append('n')
            # Select the velocity model, use LABasin
            opts.append('LABasin500')
            # Select method
            opts.append(Methods.options[method])
            # Source file
            opts.append('1')
            if method == Methods.ucsb:
                opts.append('northridge_eq_ucsb.src')
            elif method == Methods.song:
                opts.append('northridge_eq_song.src')
            else:
                opts.append('northridge_eq_gp.src')
            if method != Methods.exsim and method != Methods.csm:
                # Use rupture generator
                opts.append('y')
            # Select station from run directory
            opts.append('1')
            opts.append("%d" %
                        (stafiles.index("northridge_3_sta.stl") + 1))
            if method == Methods.exsim:
                # No custom template for ExSIM
                opts.append('n')
            if method != Methods.csm:
                # Run site response (CSM doesn't ask this question)
                opts.append('y')
            # No plots
            opts.append('n')
            opts.append('n')
            print("\t %s" % (str(opts)))
            numuser = numuser + 1
            optfiles[optfile] = opts

        print("Number of validation events: %d" % (numvalid))
        print("Number of user events: %d" % (numuser))
        print("Total aceptance tests: %d" % (numuser + numvalid))
        return optfiles

    def generate_xml(self, optfiles):
        install = InstallCfg.getInstance()

        # Generate xml workflows
        tests = []
        for key in optfiles.keys():
            sim_id = int(seqnum.get_seq_num())
            test = key.split('.')[0]
            xmlfile = os.path.join(self.input_dir, "%s.xml" % (test))
            if os.path.basename(xmlfile) in self.resume_list:
                # Skip this test
                print("Skipping %s" % (key))
                continue

            print("Generating %s" % (key))
            optfile = os.path.join(self.input_dir, key)

            # Save the option file
            op = open(optfile, 'w')
            for line in optfiles[key]:
                op.write("%s\n" % (line))
            op.close()

            # Generate xml
            print("Generating xml for %s" % (key))
            print("\t %s" % (str(optfiles[key])))
            cmd = ("%s --expert -s %d -g -o %s" %
                   (os.path.join(install.A_COMP_DIR, "run_bbp.py"),
                   sim_id, optfile))
            print("Running: %s" % (cmd))
            rc = bband_utils.runprog(cmd, False)
            if rc != 0:
                print("Failed to run bbp, aborting.")
                return []
            oldxmlfile = os.path.join(install.A_XML_DIR, "%d.xml" % (sim_id))
            shutil.copy2(oldxmlfile, xmlfile)
            if not os.path.exists(xmlfile):
                print("Workflow %s not found, aborting." % (xmlfile))
                return []

            tests.append([sim_id, xmlfile])
            time.sleep(1)

        return tests

    def run_tests(self, tests):
        install = InstallCfg.getInstance()

        # Run the tests and save results as reference data
        for test in tests:
            if os.path.basename(test[1]) in self.resume_list:
                # Skip this test
                print("Skipping %s" % (os.path.basename(test[1])))
                continue

            # Execute each test
            cmd = ("%s -s %d -x %s" %
                   (os.path.join(install.A_COMP_DIR, "run_bbp.py"),
                   test[0], test[1]))
            rc = bband_utils.runprog(cmd, False)
            if rc != 0:
                print("Failed to run acceptance test %d-%s, aborting." %
                      (test[0], test[1]))
                return 1

            # Save the bbp and rd50 files
            test_name = os.path.basename(test[1]).split('.')[0]
            cmd = "mkdir -p %s" % (os.path.join(self.ref_dir, test_name))
            bband_utils.runprog(cmd)
            rd50files = glob.glob("%s/%d/%d.*.rd50" %
                                  (install.A_OUT_DATA_DIR, test[0], test[0]))
            if len(rd50files) < 1:
                print("Did not find expected RotD50 files")
                return 1
            for rd50_file in rd50files:
                filecomps = os.path.basename(rd50_file).split('.')
                shutil.copy2(rd50_file,
                             os.path.join(self.ref_dir, test_name,
                                          "%s.rd50" % (filecomps[1])))

            # Write progress to checkpoint file
            resume_fp = open(os.path.join(install.A_OUT_LOG_DIR,
                                          "gen_resume.txt"), 'a')
            resume_fp.write("%s\n" % os.path.basename(test[1]))
            resume_fp.flush()
            resume_fp.close()

        return 0

if __name__ == '__main__':
    # Make sure BBP_DATA_DIR is set, otherwise, the produced XML files
    # will not be useable by other people
    if not 'BBP_DATA_DIR' in os.environ:
        print("Please set BBP_DATA_DIR and try again!")
        sys.exit(1)
    parser = optparse.OptionParser()
    parser.add_option("-o", "--overwrite", action="store_true",
                      dest="overwrite",
                      help="Overwrite reference solution")
    (options, args) = parser.parse_args()

    # Generate options
    generator = GenAcceptTests()
    option_list = generator.generate_option_list()
    if len(option_list) == 0:
        print("No options were produced")
        sys.exit(1)

    # Generate XML test files
    test_list = generator.generate_xml(option_list)
    if len(test_list) == 0:
        print("No tests were produced")
        sys.exit(1)

    if options.overwrite == True:
        # Run and save the solutions
        if generator.run_tests(test_list) != 0:
            print("Failed to execute acceptance tests")
            sys.exit(1)

    sys.exit(0)
