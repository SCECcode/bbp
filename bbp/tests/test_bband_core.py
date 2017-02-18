#! /usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are unit tests for the broadband platform core
$Id: test_bband_core.py 1734 2016-09-13 17:38:17Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import unittest

# Import Broadband modules
from test_bband_utils import TestBBandUtils
from test_python_code import TestPythonCode
from test_arias import TestArias
from test_bbp_format import TestBBPFormat

class CoreTestSuite(unittest.TestSuite):
    """
    Test suite containing generic platform tests
    """
    def __init__(self):
        """
        Add core Broadband functionality
        """
        unittest.TestSuite.__init__(self)
        self.addTest(unittest.makeSuite(TestBBandUtils))
        self.addTest(unittest.makeSuite(TestPythonCode))
        self.addTest(unittest.makeSuite(TestArias))
        self.addTest(unittest.makeSuite(TestBBPFormat))

if __name__ == '__main__':
    CORE = CoreTestSuite()
    unittest.TextTestRunner(verbosity=2).run(CORE)
