#! /usr/bin/python
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
