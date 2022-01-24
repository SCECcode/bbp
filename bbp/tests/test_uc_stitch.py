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
from uc_stitch import UCStitch

class TestUCStitch(unittest.TestCase):
    """
    Acceptance Test for uc_stitch.py
    """

    def setUp(self):
        self.install = InstallCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.velocity = "nr02-vs500_lf.vel"
        self.srffile = "test_ucsb.srf"
        self.stations = "test_stat.txt"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.splitext(self.stations)[0]
        self.a_tmpdir_mod = os.path.join(self.install.A_TMP_DATA_DIR,
                                         str(self.sim_id),
                                         "uc_stitch_%s" % (sta_base))

        # Create directories
        a_refdir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        cmd = "mkdir -p %s" % (a_indir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_tmpdir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_outdir)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s" % (a_logdir)
        bband_utils.runprog(cmd)

        cmd = "cp %s %s" % (os.path.join(a_refdir, self.velocity), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.srffile), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "gp", self.stations), a_indir)
        bband_utils.runprog(cmd)
        for i in range(1, 6):
            cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                             "gp", "s%02d-lf.bbp" % (i)),
                                os.path.join(a_tmpdir,
                                             "%d.s%02d-lf.bbp" %
                                             (self.sim_id, i)))
            bband_utils.runprog(cmd)
            cmd = "cp %s %s" % (os.path.join(a_refdir, "s%02d.3comp" % (i)),
                                os.path.join(a_tmpdir, "%d.s%02d.bbp" % (self.sim_id, i)))
            bband_utils.runprog(cmd)

    def test_ucsb_stitch(self):
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
                                   "%d.s%02d-stitch.bbp" % (self.sim_id, i))
            self.failIf(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                        "output merged BBP file "
                        "%s does not match reference merged bbp file %s" %
                        (bbpfile, ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestUCStitch)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
