#! /usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the broadband platforms
$Id: test_ucsb.py 1780 2017-01-10 18:17:06Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import unittest

# Import Broadband modules
from test_uc_fault_utils import TestUCFaultUtils
from test_ucrmg import TestUCrmg
from test_syn1d import TestSyn1D
from test_uc_site import TestUCSite
from test_rotd50 import TestRotD50
from test_gp_gof import TestGPGof

class UCSBTestSuite(unittest.TestSuite):
    """
    Tests for the UCSB method
    """
    def __init__(self):
        """
        Add all UCSB tests to the test suite
        """
        unittest.TestSuite.__init__(self)
        self.addTest(unittest.makeSuite(TestUCFaultUtils))
        self.addTest(unittest.makeSuite(TestUCrmg))
        self.addTest(unittest.makeSuite(TestSyn1D))
        self.addTest(unittest.makeSuite(TestUCSite))
        self.addTest(unittest.makeSuite(TestRotD50))
        self.addTest(unittest.makeSuite(TestGPGof))

if __name__ == '__main__':
    UTS = UCSBTestSuite()
    unittest.TextTestRunner(verbosity=2).run(UTS)
