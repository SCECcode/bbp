#! /usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

These are unit tests for the broadband platform core
"""
from __future__ import division, print_function

# Import Python modules
import sys
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
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(CORE)
    sys.exit(not RETURN_CODE.wasSuccessful())
