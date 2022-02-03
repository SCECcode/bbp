#!/usr/bin/env python
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
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import seqnum
import bband_utils
import arias_duration
from install_cfg import InstallCfg

class TestArias(unittest.TestCase):
    """
    Unit test for the arias duration module
    """

    def setUp(self):
        self.install = InstallCfg()
        self.sim_id = int(seqnum.get_seq_num())
        self.a_outdir = os.path.join(self.install.A_OUT_DATA_DIR,
                                     str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([self.a_outdir])

    def test_arias_duration(self):
        """
        Run the arias intensity unit test
        """
        in_file = os.path.join(self.install.A_TEST_REF_DIR, "arias",
                               "inputs", "NGA_no_1063_RRS228.AT2")
        out_file = os.path.join(self.a_outdir, "NGA_no_1063_RRS228.AT2.AD.bbp")
        ref_file = os.path.join(self.install.A_TEST_REF_DIR, "arias",
                                "reference",
                                "NGA_no_1063_RRS228.AT2.AD.bbp.ref")

        if arias_duration.ad_from_acc(in_file, out_file) != 0:
            print("Error converting ACC to Arias Duration")
            self.assertTrue(False)
        else:
            res_file = open(out_file, 'r')
            lines = res_file.readlines()
            res_file.close()

            ref_file = open(ref_file, 'r')
            rlines = ref_file.readlines()
            ref_file.close()

            # Only check if we have the same number of lines
            self.assertTrue(len(lines) == len(rlines))

if __name__ == "__main__":
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestArias)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
