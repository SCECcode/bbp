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

Unit tests for the BBplatform codes
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import bband_utils

class TestBBandUtils(unittest.TestCase):
    """
    Acceptance Test for bband_utils
    """

    def setUp(self):
        os.chdir('..')

    def tearDown(self):
        os.chdir('./tests')

    def test_runprog(self):
        """
        Test runprog function
        """
        cmd = "hostname >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        # print "returned %d" % (res)
        if not res == 0:
            self.fail("hostname command did not return success returncode")

    def test_runprog2(self):
        """
        Another test for the runprog function
        """
        cmd = "issuing_intentionally_nonexistant_command >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        if res == 0:
            self.fail("Did not properly detect an error running a command on the cmdline")

    def test_runprog3(self):
        """
        Yet another test for the runprog function
        """
        cmd = "produce >/dev/null 2>&1"
        res = bband_utils.runprog(cmd, False)
        if res == 0:
            self.fail("Sub process should returned non-zero and it didn't.")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestBBandUtils)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
