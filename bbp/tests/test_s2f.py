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

These are unit tests for the stas2files program
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

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
        self.assertFalse(filecmp.cmp(ucref, ofile) == False, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(Test_s2f)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
