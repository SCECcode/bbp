#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the SONG RMG module
$Id: test_rmg.py 1734 2016-09-13 17:38:17Z fsilva $
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
from rmg_cfg import RMGCfg
from rmg import RMG

class TestRMG(unittest.TestCase):
    """
    Unit Test for rmg.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.cfg = RMGCfg()
        os.chdir(self.install.A_COMP_DIR)
        self.sim_id = int(seqnum.get_seq_num())
        self.srcfile = "test_whittier_song.src"
        self.outsrf = "%d_test_eq.srf" % self.sim_id

        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir],
                           print_cmd=False)

        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "song",
                                         self.srcfile),
                            os.path.join(self.install.A_IN_DATA_DIR,
                                         str(self.sim_id)))
        bband_utils.runprog(cmd)

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_rmg(self):
        a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "song")
        a_res_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))

        rmg = RMG(None, self.srcfile, self.outsrf, None, sim_id=self.sim_id)
        rmg.run()
        #
        # Test conversion from RMG to srf file
        #

        a_ref_file = os.path.join(a_ref_dir, "song_test_whittier.srf")
        a_newfile = os.path.join(a_res_dir, self.outsrf)
        errmsg = ("Output file %s does not match reference file %s" %
                  (a_newfile, a_ref_file))
        self.failIf(not cmp_bbp.cmp_srf(a_ref_file, a_newfile,
                                        tolerance=0.0011) == 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestRMG)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
