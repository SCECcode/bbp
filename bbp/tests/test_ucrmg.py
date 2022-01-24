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
"""
from __future__ import division, print_function

# Import Python modules
import os
import shutil
import unittest

# Import Broadband modules
import seqnum
import bband_utils
import cmp_bbp
from install_cfg import InstallCfg
from ucrmg_cfg import UCrmgCfg
from ucrmg import UCrmg

class TestUCrmg(unittest.TestCase):
    """
    Acceptance Test for URrmg. This assumes the acceptance test calling
    script runs the basic executable.
    """

    def setUp(self):
        self.install = InstallCfg()
        self.r_velmodel = "nr02-vs500_lf.vel"
        self.vmodel_name = "LABasin500"
        self.r_srcfile = "test_wh_ucsb.src"
        self.r_srffile = "FFSP_OUTPUT.bst"
        self.sim_id = int(seqnum.get_seq_num())
        self.cfg = UCrmgCfg(self.vmodel_name)

        a_refdir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir, a_logdir],
                           print_cmd=False)

        # Copy files
        shutil.copy2(os.path.join(a_refdir, self.r_srcfile),
                     os.path.join(a_indir, self.r_srcfile))
        shutil.copy2(os.path.join(a_refdir, self.r_velmodel),
                     os.path.join(a_indir, self.r_velmodel))

        os.chdir(a_tmpdir)

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_ucgen(self):
        """
        Unit test for the UCSB rupture generator
        """
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_res_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))

        a_ref_file = os.path.join(a_ref_dir, self.r_srffile)
        a_newfile = os.path.join(a_res_dir, self.r_srffile)

        uc_obj = UCrmg(self.r_velmodel, self.r_srcfile,
                       self.r_srffile, self.vmodel_name,
                       sim_id=self.sim_id)
        uc_obj.run()

        errmsg = "Output file does not match reference file: %s" % (a_newfile)
        self.assertFalse(not cmp_bbp.cmp_ffsp(a_ref_file,
                                              a_newfile,
                                              tolerance=0.0025) == 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestUCrmg)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
