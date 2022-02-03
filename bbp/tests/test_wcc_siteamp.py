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
from wcc_siteamp_cfg import WccSiteampCfg
from wcc_siteamp import WccSiteamp

class TestWccSiteamp(unittest.TestCase):
    """
    Acceptance Test for wcc_siteamp.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.wcc_siteamp_cfg = WccSiteampCfg("LABasin500", "GP")
        os.chdir(self.install.A_INSTALL_ROOT)
        self.stations = "test_stat.txt"
        self.sim_id = int(seqnum.get_seq_num())
        self.freqs = ['lf', 'hf']

        refdir = os.path.join(self.install.A_TEST_REF_DIR, "gp")
        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir], print_cmd=False)

        # Copy station list
        cmd = "cp %s %s" % (os.path.join(refdir, self.stations), indir)
        bband_utils.runprog(cmd)

        # Copy LF and HF seismograms
        for i in range(1, 6):
            for freq in self.freqs:
                cmd = "cp %s %s" % (os.path.join(refdir,
                                                 "s%02d-%s.bbp" %
                                                 (i, freq)),
                                    os.path.join(tmpdir,
                                                 "%d.s%02d-%s.bbp" %
                                                 (self.sim_id, i, freq)))
                bband_utils.runprog(cmd)

    def test_wcc_siteamp(self):
        """
        Test GP site response module
        """
        wcc_obj = WccSiteamp(self.stations, "GP",
                             "LABasin500", sim_id=self.sim_id)
        wcc_obj.run()
        freqs = ['lf', 'hf']

        for i in range(1, 6):
            for freq in freqs:
                ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                        "gp",
                                        "s%02d-%s-site.bbp" %
                                        (i, freq))
                bbpfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                       str(self.sim_id),
                                       "%d.s%02d-%s.acc.bbp" %
                                       (self.sim_id, i, freq))
                self.assertFalse(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                                 "output %s BBP file %s does not match reference %s bbp file %s " %
                                 (freq, bbpfile, freq, ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestWccSiteamp)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
