#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This is the unit test for the as16.py BBP module
$Id: test_as16.py 1795 2017-02-09 16:23:34Z fsilva $
"""

# Import Python modules
import os
import sys
import unittest

# Import Broadband modules
import cmp_bbp
import seqnum
import bband_utils
from install_cfg import InstallCfg
from as16 import AS16
from station_list import StationList

class TestAS16(unittest.TestCase):
    """
    Unit test for as16.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.stations = "nr_v13_3_1.stl"
        self.source = "nr_v14_02_1.src"        
        self.eventname = "NR"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])
        sim_id = self.sim_id

        # Set up paths
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(self.install.A_TMP_DATA_DIR,
                                     str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir,
                                           "validations",
                                           "stewart_duration_gmpe")

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_logdir, a_validation_outdir])

        # Copy station list
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "as16", self.stations),
                            a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "as16", self.source),
                            a_indir)
        bband_utils.runprog(cmd)

    def test_as16(self):
        """
        Run the AS16 GMPE test
        """
        as16_obj = AS16(self.stations, self.source,
                        self.eventname, sim_id=self.sim_id)
        as16_obj.run()

        # Check results
        ref_sum_file = os.path.join(self.install.A_TEST_REF_DIR, "as16",
                                    "as16.%s.txt" % (self.eventname))
        cal_sum_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                    str(self.sim_id),
                                    "validations",
                                    "stewart_duration_gmpe",
                                    "%d.as16.%s.txt" %
                                    (self.sim_id, self.eventname))
        self.failIf(cmp_bbp.cmp_anderson_gof(ref_sum_file, cal_sum_file,
                                             tolerance=0.005,
                                             start_col=1,
                                             sep=",") != 0,
                    "AS16 Summary file does not match reference file!")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sim_id = int(sys.argv[1])
        print "Will test with sim_id: %d" % (sim_id)
    else:
        sim_id = None
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAS16)
    unittest.TextTestRunner(verbosity=2).run(suite)
