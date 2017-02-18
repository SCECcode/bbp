#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for syn1D_LAH.py
$Id: test_syn1d.py 1734 2016-09-13 17:38:17Z fsilva $
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
from syn1D_cfg import Syn1DCfg
from syn1D import Syn1D

class TestSyn1D(unittest.TestCase):
    """
    Acceptance Test for syn1D. This assumes the acceptance test calling
    script runs the basic executable.
    """

    def setUp(self):
        self.r_velmodel = "labasin.vel"
        self.r_srcfile = "test_wh_ucsb.src"
        self.r_metadata = "metadata.txt"
        self.r_stations = "one_stat.txt"
        self.r_srffile = "FFSP_OUTPUT.001"
        self.vmodel_name = "LABasin"
        self.sim_id = int(seqnum.get_seq_num())

        self.install = InstallCfg()
        self.cfg = Syn1DCfg(self.vmodel_name)

        cmd = "mkdir -p %s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_LOG_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.r_srffile,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
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
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.r_stations,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
        bband_utils.runprog(cmd)
        cmd = ("cp %s/ucsb/faultGlobal.in %s/%d/." %
               (self.install.A_TEST_REF_DIR, self.install.A_TMP_DATA_DIR,
                self.sim_id))
        bband_utils.runprog(cmd)
        cmd = ("cp %s/ucsb/source_model.list %s/%d/." %
               (self.install.A_TEST_REF_DIR, self.install.A_TMP_DATA_DIR,
                self.sim_id))
        bband_utils.runprog(cmd)

        os.chdir(os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id)))

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_syn1d(self):
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_res_dir = os.path.join(self.install.A_TMP_DATA_DIR ,str(self.sim_id))

        syn_obj = Syn1D(self.r_velmodel, self.r_srcfile, self.r_srffile,
                        self.r_stations, self.vmodel_name, sim_id=self.sim_id)
        syn_obj.run()

        # Check each station
        for i in range(1, 2):
            a_ref_file = os.path.join(a_ref_dir,
                                      "s%02d.3comp" % (i))
            a_newfile = os.path.join(a_res_dir,
                                     "stitch_one_stat",
                                     "s%02d.3comp" % (i))

            errmsg = ("Output file %s does not match reference file: %s" %
                      (a_newfile, a_ref_file))
            self.failIf(cmp_bbp.cmp_bbp(a_ref_file, a_newfile) != 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestSyn1D)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
