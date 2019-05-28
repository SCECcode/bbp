#! /usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

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
import cmp_bbp
import seqnum
import bband_utils
from install_cfg import InstallCfg
from syn1D import Syn1D

class TestSyn1D(unittest.TestCase):
    """
    Acceptance Test for syn1D. This assumes the acceptance test calling
    script runs the basic executable.
    """

    def setUp(self):
        self.r_velmodel = "nr02-vs500_lf.vel"
        self.r_srcfile = "test_wh_ucsb.src"
        self.r_stations = "one_stat.txt"
        self.r_srffile = "FFSP_OUTPUT.001"
        self.vmodel_name = "LABasin500"
        self.sim_id = int(seqnum.get_seq_num())

        self.install = InstallCfg()

        a_refdir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))

        cmd = "mkdir -p %s" % (a_indir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_tmpdir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_outdir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_logdir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.r_srffile), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.r_srcfile), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.r_velmodel), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.r_stations), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, "source_model.list"), a_tmpdir)
        bband_utils.runprog(cmd)

        os.chdir(a_tmpdir)

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_syn1d(self):
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_res_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))

        syn_obj = Syn1D(self.r_velmodel, self.r_srcfile, self.r_srffile,
                        self.r_stations, self.vmodel_name, sim_id=self.sim_id)
        syn_obj.run()

        # Check each station
        for i in range(1, 2):
            a_ref_file = os.path.join(a_ref_dir,
                                      "s%02d.3comp" % (i))
            a_newfile = os.path.join(a_res_dir,
                                     "stitch_%s" %
                                     (os.path.splitext(self.r_stations)[0]),
                                     "s%02d.3comp" % (i))

            errmsg = ("Output file %s does not match reference file: %s" %
                      (a_newfile, a_ref_file))
            self.failIf(cmp_bbp.cmp_bbp(a_ref_file, a_newfile) != 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestSyn1D)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
