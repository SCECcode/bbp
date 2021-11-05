#! /usr/bin/env python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

These are unit test for the Broadband Python code
"""
from __future__ import division, print_function

# Import Python modules
import os
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
        cmd = ("python %s -v >/dev/null" %
               (os.path.join(self.install.A_COMP_DIR,
                             "run_bbp.py")))
        self.assertFalse(bband_utils.runprog(cmd, False) != 0,
                         "Cannot start Broadband plotform!")

    def test_python_code_comps(self):
        """
        Run Python with -tt flag to detect mix of tabs and spaces in the code
        """
        self.install = InstallCfg()
        cmd = ("python -tt -m compileall -f -q -l %s" %
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
    unittest.TextTestRunner(verbosity=2).run(SUITE)
