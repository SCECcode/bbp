#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the jbsim.py
$Id: test_uc_fault_utils.py 1734 2016-09-13 17:38:17Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import filecmp
import unittest

# Import Broadband modules
import seqnum
import bband_utils

from install_cfg import InstallCfg
import uc_fault_utils

class TestUCFaultUtils(unittest.TestCase):
    """
    Acceptance Test for jbsim.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.r_srcfile = "test_wh_ucsb.src"
        self.r_faultfile = "ffsp.inp"
        self.r_velmodel = "labasin.vel"
        self.vmodel_name = "LABasin"
        self.sim_id = int(seqnum.get_seq_num())
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

    def test_uc_fault_utils(self):
        """
        Test UCSB fault utilities
        """
        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_src_file = os.path.join(indir, self.r_srcfile)
        a_ffsp_inp = os.path.join(tmpdir, self.r_faultfile)

        # Create ffsp.inp file
        uc_fault_utils.uc_create_ffsp_inp(a_ffsp_inp, a_src_file,
                                          self.vmodel_name)
        ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                "ucsb", self.r_faultfile)
        self.failIf(filecmp.cmp(ref_file, a_ffsp_inp) == False,
                    "output fault file does not match reference fault file")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestUCFaultUtils)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
