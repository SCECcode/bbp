#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are unit test for the Broadband Python code
$Id: test_python_code.py 1734 2016-09-13 17:38:17Z fsilva $
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
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Cannot start Broadband plotform!")

    def test_python_code_comps(self):
        """
        Run Python with -tt flag to detect mix of tabs and spaces in the code
        """
        self.install = InstallCfg()
        cmd = ("python -tt -m compileall -f -q -l %s" %
               (self.install.A_COMP_DIR))
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Python code in comps directory mixes tabs and spaces!")

    def test_python_code_tests(self):
        """
        Run Python with -tt flag to detect mix of tabs and spaces in the code
        """
        self.install = InstallCfg()
        cmd = ("python -tt -m compileall -f -q -l %s" %
               (self.install.A_TEST_DIR))
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Python code in test directory mixes tabs and spaces!")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestPythonCode)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
