#! /usr/bin/python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

These are unit tests for the coordiante conversion codes
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
        self.assertFalse(filecmp.cmp(reffile, ofile) == False,
                         "LL to XY did not work")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestCC)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
