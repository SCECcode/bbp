#!/usr/bin/env python
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

These are acceptance tests for the broadband platforms
"""
from __future__ import division, print_function

# Import Python modules
import os
import new
import sys
import shutil
import optparse
import unittest

# Import Broadband modules
import bband_utils
import seqnum
import cmp_bbp
from install_cfg import InstallCfg

def find_tests(test, rerun):
    """
    # This function searches for .xml files in the accept_inputs directory
    """

    install = InstallCfg()
    resume = True

    accept_test_inputs = "accept_inputs"
    accept_test_refs = "accept_refs"

    input_dir = os.path.join(install.A_TEST_REF_DIR, accept_test_inputs)
    if not os.path.exists(input_dir):
        # These are expected to be in the dist
        print("Acceptance test inputs dir %s does not exist, aborting" %
              (input_dir))
        sys.exit()

    # Create list of test XML files
    files = os.listdir(input_dir)
    wfext = ".xml"

    # First we find all the tests
    test_files = []
    for testfile in files:
        if testfile.endswith(wfext):
            # Don't add SDSU tests on Mac OS X
            if sys.platform == 'darwin' and testfile.find("SDSU") >= 0:
                if test is None or (test is not None and testfile.find(test) >= 0):
                    print("*** Mac OS X detected: skipping test %s." %
                          (testfile))
                continue
            if test is None:
                test_files.append(testfile)
            else:
                if testfile.find(test) >= 0:
                    test_files.append(testfile)

    resume_file = os.path.join(install.A_OUT_LOG_DIR, "resume.txt")
    resume_list = ""
    if rerun:
        os.remove(resume_file)
    # Check for already completed tests if not rerunning
    if resume == True and rerun == False:
        if os.path.exists(resume_file):
            resume_fp = open(resume_file, 'r')
            resume_list = resume_fp.read().splitlines()
            completed_test_count = len(resume_list)
            print("==> Completed Tests : %d" % (completed_test_count))
            resume_fp.close()
            if ((test is None) and
                (completed_test_count >= len(test_files))):
                print("All the acceptance tests have passed previously!")
                proceed = raw_input("Would you like to re-run "
                                    "all the acceptance tests? (y/n)")
                if str.lower(proceed) == 'y':
                    os.remove(resume_file)
                    resume_list = ""
                else:
                    sys.exit(0)

    # Create unittest test case for each file
    for xml_file in test_files:

        # Skip test if we ran it already
        if xml_file in resume_list:
            print("==> Skipping %s" % (xml_file))
            continue

        file_base = xml_file[0:xml_file.find(wfext)]
        # pieces = file_base.split('-')

        # Adjust tolerance depending on test mode
        tolerance = 0.03

        #This defines a method that we're going to add to the
        #BBPAcceptanceTests class. The keyword binding has to
        #be done b/c Python is storing pointers to 'file' and 'file_base'
        #so w/o the keywords, 'file' and 'file_base' in the function will
        #point to the final values
        def permutation_test(self, file_base=file_base, xml_file=xml_file):
            input_dir = os.path.join(self.install.A_TEST_REF_DIR,
                                     accept_test_inputs)
            log_dir = os.path.join(self.install.A_OUT_LOG_DIR,
                                   "acceptance_test_logs")
            sim_id = int(seqnum.get_seq_num())
            self.file_base = file_base
            self.log_file = os.path.join(log_dir, "%s.log" % (self.file_base))
            self.input_file = os.path.join(input_dir, xml_file)
            cmd = ("%s/run_bbp.py -x %s -s %d -l %s" %
                   (self.install.A_COMP_DIR,
                    self.input_file, sim_id, self.log_file))
            rc = bband_utils.runprog(cmd, False)
            self.failIf(rc != 0, "Acceptance test failed to execute")
            ref_file_dir = os.path.join(self.install.A_TEST_REF_DIR,
                                        accept_test_refs,
                                        self.file_base)
            agree = True
            for ref_file in os.listdir(ref_file_dir):
                if os.path.isfile(os.path.join(ref_file_dir, ref_file)):
                    test_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                             str(sim_id),
                                             ("%d.%s" % (sim_id, ref_file)))
                    a_ref_file = os.path.join(ref_file_dir, ref_file)
                    compare_result = cmp_bbp.cmp_bbp(a_ref_file, test_file,
                                                     tolerance=tolerance)
                    errmsg = ("Output file "
                              "%s does not match reference file: %s" %
                              (test_file, a_ref_file))
                    self.failIf(compare_result != 0, errmsg)
                    if compare_result != 0:
                        agree = False
            if agree == True:
                # Write success to the resume file
                resume_fp = open(os.path.join(install.A_OUT_LOG_DIR,
                                              "resume.txt"), 'a')
                resume_fp.write("%s\n" % xml_file)
                resume_fp.flush()
                resume_fp.close()
            sys.stdout.flush()
            sys.stderr.flush()

        # We create a method object which is an instance method for
        # BBPAcceptanceTests which executes the code in
        # testPermutation
        method = new.instancemethod(permutation_test,
                                    None, BBPAcceptanceTests)
        # We give the method a new name in BBPAcceptanceTests
        # which contains the xml file being run
        setattr(BBPAcceptanceTests, "test_%s" % file_base, method)

class BBPAcceptanceTests(unittest.TestCase):

    def setUp(self):
        self.install = InstallCfg()
        accept_test_inputs = "accept_inputs"
        src_path = ""
        self.resume = True
        run_dir = self.install.A_USER_DATA_DIR

        # Create run directory, in case it doesn't exist
        bband_utils.mkdirs([run_dir], print_cmd=False)

        if not os.path.exists(os.path.join(run_dir, "northridge_3_sta.stl")):
            src_path = os.path.join(self.install.A_TEST_REF_DIR,
                                    accept_test_inputs,
                                    "northridge_3_sta.stl")
            shutil.copy2(src_path, run_dir)

        if not os.path.exists(os.path.join(run_dir, "northridge_eq_gp.src")):
            src_path = os.path.join(self.install.A_TEST_REF_DIR,
                                    accept_test_inputs,
                                    "northridge_eq_gp.src")
            shutil.copy2(src_path, run_dir)

        if not os.path.exists(os.path.join(run_dir, "northridge_eq_ucsb.src")):
            src_path = os.path.join(self.install.A_TEST_REF_DIR,
                                    accept_test_inputs,
                                    "northridge_eq_ucsb.src")
            shutil.copy2(src_path, run_dir)

        if not os.path.exists(os.path.join(run_dir, "northridge_eq_song.src")):
            src_path = os.path.join(self.install.A_TEST_REF_DIR,
                                    accept_test_inputs,
                                    "northridge_eq_song.src")
            shutil.copy2(src_path, run_dir)

        if not os.path.exists(os.path.join(self.install.A_OUT_LOG_DIR,
                                           "acceptance_test_logs")):
            bband_utils.mkdirs([os.path.join(self.install.A_OUT_LOG_DIR,
                                             "acceptance_test_logs")])

if __name__ == '__main__':
    # Parse options
    parser = optparse.OptionParser()
    parser.add_option("-t", "--test",
                      dest="test",
                      help="Execute specific test",
                      metavar="TEST")
    parser.add_option("-r", "--rerun",
                      action="store_true",
                      dest="rerun",
                      help="Rerun tests already completed")

    (options, args) = parser.parse_args()

    if options.test is not None:
        test = options.test
    else:
        test = None
    if options.rerun is not None:
        rerun = True
    else:
        rerun = False

    find_tests(test, rerun)
    suite = unittest.TestLoader().loadTestsFromTestCase(BBPAcceptanceTests)
    print("==> Number of tests to run: %d" % suite.countTestCases())
    unittest.TextTestRunner(verbosity=2).run(suite)
