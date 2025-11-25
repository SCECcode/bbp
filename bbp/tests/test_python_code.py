#! /usr/bin/env python3
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

These are unit test for the Broadband Python code
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import unittest

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg

class TestPythonCode(unittest.TestCase):
    """
    Unit test for Broadband's Python code
    """
    def test_execute_platform_bbp(self):
        """
        Run Broadband Plotform to make sure we can start it
        """
        self.install = InstallCfg()
        cmd = ("python3 %s -v >/dev/null" %
               (os.path.join(self.install.A_COMP_DIR,
                             "run_bbp.py")))
        self.assertFalse(bband_utils.runprog(cmd, False) != 0,
                         "Cannot start Broadband plotform!")

    def test_python_code_comps(self):
        """
        Run Python with -tt flag to detect mix of tabs and spaces in the code
        """
        self.install = InstallCfg()
        cmd = ("python3 -tt -m compileall -f -q -l %s" %
               (self.install.A_COMP_DIR))
        self.assertFalse(bband_utils.runprog(cmd, False) != 0,
                         "Python code in comps directory mixes tabs and spaces!")

    def test_python_code_tests(self):
        """
        Run Python with -tt flag to detect mix of tabs and spaces in the code
        """
        self.install = InstallCfg()
        cmd = ("python -tt -m compileall -f -q -l %s" %
               (self.install.A_TEST_DIR))
        self.assertFalse(bband_utils.runprog(cmd, False) != 0,
                         "Python code in test directory mixes tabs and spaces!")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestPythonCode)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
