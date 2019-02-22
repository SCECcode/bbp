#! /usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Unit tests for the BBplatform codes
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import bband_utils

class TestBBandUtils(unittest.TestCase):
    """
    Acceptance Test for bband_utils
    """

    def setUp(self):
        os.chdir('..')

    def tearDown(self):
        os.chdir('./tests')

    def test_runprog(self):
        """
        Test runprog function
        """
        cmd = "hostname >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        # print "returned %d" % (res)
        if not res == 0:
            self.fail("hostname command did not return success returncode")

    def test_runprog2(self):
        """
        Another test for the runprog function
        """
        cmd = "issuing_intentionally_nonexistant_command >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        if res == 0:
            self.fail("Did not properly detect an error running a command on the cmdline")

    def test_runprog3(self):
        """
        Yet another test for the runprog function
        """
        cmd = "produce >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        if res == 0:
            self.fail("Sub process should returned non-zero and it didn't.")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestBBandUtils)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
