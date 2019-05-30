#! /usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
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
