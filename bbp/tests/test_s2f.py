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

These are unit tests for the stas2files program
"""
from __future__ import division, print_function

# Import Python modules
import os

# Import BBP modules
import seqnum
import filecmp
import unittest
import stas2files
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList

class Test_s2f(unittest.TestCase):
    """
    This tests the station format conversion code
    """
    def test_gp2uc(self):
        """
        Inputs a GP format file and get out a UCSB format file
        """
        install = InstallCfg()
        sim_id = int(seqnum.get_seq_num())
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_refdir = os.path.join(install.A_TEST_REF_DIR, "ucsb")

        #
        # Make sure output directories exist
        #
        bband_utils.mkdirs([a_tmpdir], print_cmd=False)

        # File paths
        gpfile = os.path.join(a_refdir, "stats-h0.125.ll")
        ucref = os.path.join(a_refdir, "stations.ll")
        ofile = os.path.join(a_tmpdir, "stations.ll")
        ofile30 = os.path.join(a_tmpdir, "stations.vs30")

        sl = StationList(gpfile)
        _ = stas2files.gp2uc_stalist(sl, ofile, ofile30)
        errmsg = "Conversion of station list from GP to UC format failed"
        self.failIf(filecmp.cmp(ucref, ofile) == False, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(Test_s2f)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
