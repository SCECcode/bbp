#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for match.py
$Id: test_uc_stitch.py 1734 2016-09-13 17:38:17Z fsilva $
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
from uc_stitch import UCStitch

class TestUCStitch(unittest.TestCase):
    """
    Acceptance Test for uc_stitch.py
    """

    def setUp(self):
        self.install = InstallCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.velocity = "labasin.vel"
        self.srffile = "test_ucsb.srf"
        self.stations = "test_stat.txt"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.splitext(self.stations)[0]
        self.a_tmpdir_mod = os.path.join(self.install.A_TMP_DATA_DIR,
                                         str(self.sim_id),
                                         "uc_stitch_%s" % (sta_base))
        cmd = "mkdir -p %s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_LOG_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.velocity,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/ucsb/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                         self.srffile,
                                         self.install.A_IN_DATA_DIR,
                                         self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.stations,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd)
        for i in range(1, 6):
            cmd = ("cp %s/gp/s%02d-lf.bbp %s/%d/%d.s%02d-lf.bbp" %
                   (self.install.A_TEST_REF_DIR, i,
                    self.install.A_TMP_DATA_DIR, self.sim_id,
                    self.sim_id, i))
            bband_utils.runprog(cmd)
            cmd = ("cp %s/ucsb/s%02d.3comp %s/%d/%d.s%02d.bbp" %
                   (self.install.A_TEST_REF_DIR, i,
                    self.install.A_TMP_DATA_DIR, self.sim_id,
                    self.sim_id, i))
            bband_utils.runprog(cmd)

    def test_bbp_wid(self):
        """
        Test the UCSB stitch code
        """
        stitch_obj = UCStitch(self.velocity, "", self.srffile,
                              self.stations, acc=False,
                              sim_id=self.sim_id)
        stitch_obj.run()

        for i in range(1, 6):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "ucsb",
                                    "s%02d.stitched.bbp" % (i))
            bbpfile = os.path.join(self.a_tmpdir_mod,
                                   "%d.s%02d-stitch.bbp" % (str(self.sim_id),
                                                            i))
            self.failIf(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                        "output merged BBP file "
                        "%s does not match reference merged bbp file %s" %
                        (bbpfile, ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestUCStitch)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
