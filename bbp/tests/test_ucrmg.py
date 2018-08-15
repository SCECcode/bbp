#! /usr/bin/env python
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
        self.r_velmodel = "labasin.vel"
        self.vmodel_name = "LABasin863"
        self.r_srcfile = "test_wh_ucsb.src"
        self.r_srffile = "FFSP_OUTPUT.001"
        self.sim_id = int(seqnum.get_seq_num())
        self.cfg = UCrmgCfg(self.vmodel_name)
        cmd = "mkdir -p %s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_LOG_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.r_srcfile,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.r_velmodel,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
        bband_utils.runprog(cmd)
        os.chdir(os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id)))

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
        self.failIf(not cmp_bbp.cmp_ffsp(a_ref_file,
                                         a_newfile,
                                         tolerance=0.0025) == 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestUCrmg)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
