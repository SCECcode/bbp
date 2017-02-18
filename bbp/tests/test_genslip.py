#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the bbcoda2.py
$Id: test_genslip.py 1740 2016-09-13 19:25:14Z fsilva $
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
from genslip_cfg import GenslipCfg
from genslip import Genslip

class TestGenslip(unittest.TestCase):
    """
    Acceptance Test for genslip.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.cfg = GenslipCfg()
        os.chdir(self.install.A_COMP_DIR)
        self.sim_id = int(seqnum.get_seq_num())
        self.velmodel = "genslip_nr_generic1d-gp01.vmod"
        self.srcfile = "test_wh.src"
        self.outsrf = "%d_test_eq.srf" % self.sim_id

        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir],
                           print_cmd=False)

        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.velmodel,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.srcfile,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd, print_cmd=False)

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_gensrf(self):
        """
        Test GP rupture generator
        """
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "gp")
        a_res_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))

        gen = Genslip(self.velmodel, self.srcfile,
                      self.outsrf, "LABasin",
                      sim_id=self.sim_id)
        gen.run()
        #
        # Test conversion from genslip to srf file
        #

        a_ref_file = os.path.join(a_ref_dir,
                                  "m5.89-0.20x0.20_s2379646.srf")
        a_newfile = os.path.join(a_res_dir, self.outsrf)
        errmsg = ("Output file %s does not match reference file %s" %
                  (a_newfile, a_ref_file))
        self.failIf(not cmp_bbp.cmp_srf(a_ref_file, a_newfile,
                                        tolerance=0.0011) == 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestGenslip)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
