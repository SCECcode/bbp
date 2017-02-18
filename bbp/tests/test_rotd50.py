#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Created on Aug 10, 2012
@author: maechlin
$Id: test_rotd50.py 1734 2016-09-13 17:38:17Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import rotd50
import seqnum
import bband_utils
import install_cfg

class TestRotD50(unittest.TestCase):
    """
    Unit test for the rotd50.py module
    """

    def setUp(self):
        """
        Sets up the environment for the test
        """
        self.install = install_cfg.InstallCfg.getInstance()
        self.sim_id = int(seqnum.get_seq_num())

        # Make sure all directories exist
        self.indir = os.path.join(self.install.A_IN_DATA_DIR,
                                  str(self.sim_id))
        self.tmpdir = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id))
        self.outdir = os.path.join(self.install.A_OUT_DATA_DIR,
                                   str(self.sim_id))
        self.logdir = os.path.join(self.install.A_OUT_LOG_DIR,
                                   str(self.sim_id))
        bband_utils.mkdirs([self.indir, self.tmpdir, self.outdir, self.logdir],
                           print_cmd=False)

    def test_rotd50(self):
        """
        Test the rotd50 module
        """
        # Load configuration, set sim_id
        sim_id = self.sim_id
        # Reference directory
        ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucb")

        r_station_list = "rotd50.stl"
        a_src_station_list = os.path.join(ref_dir, r_station_list)
        a_dst_station_list = os.path.join(self.indir, r_station_list)

        # Copy station list to indir
        cmd = "cp %s %s" % (a_src_station_list, a_dst_station_list)
        bband_utils.runprog(cmd)

        # Copy velocity bbp files to outdir
        for i in range(1, 6):
            src_bbp = os.path.join(ref_dir, "8000%d.vel.bbp" % (i))
            dst_bbp = os.path.join(self.outdir, "%d.8000%d.vel.bbp" %
                                   (sim_id, i))
            cmd = "cp %s %s" % (src_bbp, dst_bbp)
            bband_utils.runprog(cmd)

        rotd50_obj = rotd50.RotD50(r_station_list, sim_id)
        rotd50_obj.run()

if __name__ == "__main__":
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestRotD50)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
