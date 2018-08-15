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
