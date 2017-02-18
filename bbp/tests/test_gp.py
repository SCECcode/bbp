#! /usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the broadband platforms
$Id: test_gp.py 1780 2017-01-10 18:17:06Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import unittest

# Import Broadband modules
from test_genslip import TestGenslip
from test_jbsim import TestJbsim
from test_hfsims import TestHfsims
from test_wcc_siteamp import TestWccSiteamp
from test_match import TestMatch
from test_respect import TestRespect
from test_rotd50 import TestRotD50
from test_gp_gof import TestGPGof

class GPTestSuite(unittest.TestSuite):
    """
    Tests for the Graves & Pitarka method
    """
    def __init__(self):
        """
        Add all GP tests to the test suite
        """
        unittest.TestSuite.__init__(self)
        self.addTest(unittest.makeSuite(TestGenslip))
        self.addTest(unittest.makeSuite(TestJbsim))
        self.addTest(unittest.makeSuite(TestHfsims))
        self.addTest(unittest.makeSuite(TestWccSiteamp))
        self.addTest(unittest.makeSuite(TestMatch))
        self.addTest(unittest.makeSuite(TestRotD50))
        self.addTest(unittest.makeSuite(TestGPGof))

if __name__ == '__main__':
    UTS = GPTestSuite()
    unittest.TextTestRunner(verbosity=2).run(UTS)
