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
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import seqnum
import bband_utils
import cmp_bbp
from install_cfg import InstallCfg
from irikura_gen_srf import IrikuraGenSrf

class TestGenSRF(unittest.TestCase):
    """
    Acceptance Test for irikura_gen_srf.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.srcfile = "whittier_v12_11_0_fs.src"
        self.outsrf = "whittier_v12_11_0_fs.srf"
        self.velmodel = "nr02-vs500.fk1d"
        self.sim_id = int(seqnum.get_seq_num())

        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        refdir = os.path.join(self.install.A_TEST_REF_DIR, "irikura")

        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir],
                           print_cmd=False)

        # Copy input files
        cmd = "cp %s %s" % (os.path.join(refdir, self.velmodel), indir)
        bband_utils.runprog(cmd, print_cmd=False)

        cmd = "cp %s %s" % (os.path.join(refdir, self.srcfile), indir)
        bband_utils.runprog(cmd, print_cmd=False)

        os.chdir(tmpdir)

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_gensrf(self):
        """
        Test Irikura rupture generator
        """
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "irikura")
        a_res_dir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))

        # Run rupture generator
        gen_srf = IrikuraGenSrf(self.velmodel, self.srcfile,
                            self.outsrf, "LABasin500",
                            sim_id=self.sim_id)
        gen_srf.run()

        #
        # Check results
        #
        a_ref_file = os.path.join(a_ref_dir, self.outsrf)
        a_newfile = os.path.join(a_res_dir, self.outsrf)
        errmsg = ("Output file %s does not match reference file %s" %
                  (a_newfile, a_ref_file))
        self.assertFalse(not cmp_bbp.cmp_srf(a_ref_file, a_newfile,
                                             tolerance=0.0011) == 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestGenSRF)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
