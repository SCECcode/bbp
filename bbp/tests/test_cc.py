#! /usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are unit tests for the coordiante conversion codes
$Id: test_cc.py 1734 2016-09-13 17:38:17Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import filecmp
import unittest

# Import Broadband modules
import cc
import seqnum
import bband_utils
from install_cfg import InstallCfg

class TestCC(unittest.TestCase):
    """
    This module tests the coordinate convertion functions in cc.py
    """
    def test_xy2ll(self):
        """
        ll2xy mlon=-118 mlat=34 xazim=0 < stats.ll > stats_out.xy
        This will return a file with a different suffix to identify the contents
        """
        self.install = InstallCfg.getInstance()
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
        bband_utils.mkdirs([self.indir, self.tmpdir, self.outdir, self.logdir])

        ilon = -118.0
        ilat = 34.0
        iaz = 0.0
        infile = os.path.join(self.indir, "100.bob.ll")
        ofile = os.path.join(self.outdir, "100.bob.xy")
        reffile = os.path.join(self.install.A_TEST_REF_DIR,
                               "sdsu", "100.bob.xy")
        # Write input file
        in_file = open(infile, "w")
        data = "%f %f\n" % (ilon, ilat)
        in_file.write(data)
        in_file.close()

        # Run the test
        cc.ll2xy(infile, ofile, ilon, ilat, iaz)

        # Check output
        self.failIf(filecmp.cmp(reffile, ofile) == False,
                    "LL to XY did not work")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestCC)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
