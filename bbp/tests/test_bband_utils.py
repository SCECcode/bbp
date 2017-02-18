#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Unit tests for the BBplatform codes
$Id: test_bband_utils.py 1734 2016-09-13 17:38:17Z fsilva $
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
