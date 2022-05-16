#! /usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

These are unit tests for the coordiante conversion codes
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
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
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
