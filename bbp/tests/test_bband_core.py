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

These are unit tests for the broadband platform core
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
