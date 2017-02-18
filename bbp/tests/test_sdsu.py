#! /usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the broadband platforms
$Id: test_sdsu.py 1780 2017-01-10 18:17:06Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import unittest

# Import Broadband modules
from test_bbtoolbox import TestBBToolbox
from test_wcc_siteamp import TestWccSiteamp
from test_rotd50 import TestRotD50
from test_gp_gof import TestGPGof
from test_sdsu_mogof import TestSDSUMOGof

class SDSUTestSuite(unittest.TestSuite):
    """
    Tests for the SDSU method
    """
    def __init__(self):
        """
        Adds individual tests to the suite
        """
        unittest.TestSuite.__init__(self)
        self.addTest(unittest.makeSuite(TestBBToolbox))
        self.addTest(unittest.makeSuite(TestWccSiteamp))
        self.addTest(unittest.makeSuite(TestRotD50))
        self.addTest(unittest.makeSuite(TestGPGof))
        self.addTest(unittest.makeSuite(TestSDSUMOGof))

if __name__ == '__main__':
    STS = SDSUTestSuite()
    unittest.TextTestRunner(verbosity=2).run(STS)
