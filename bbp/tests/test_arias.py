#!/usr/bin/env python
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
