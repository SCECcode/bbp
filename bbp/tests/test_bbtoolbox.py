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
from bbtoolbox_cfg import BBToolboxCfg
from bbtoolbox import BBToolbox

class TestBBToolbox(unittest.TestCase):
    """
    Contans unit test for SDSU BBToolbox module
    """

    def setUp(self):
        self.install = InstallCfg()
        self.cfg = BBToolboxCfg()
        self.sim_id = int(seqnum.get_seq_num())
        self.velmodel = "sdsu-aug2018-labasin-vmod.txt"
        self.srffile = "m589-s2379646.srf"
        self.stations = "test_stat.txt"
        self.srcfile = "wh_test.src"
        self.vmodel_name = "LABasin500"

        # Set up paths
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        a_refdir = os.path.join(self.install.A_TEST_REF_DIR, "sdsu")

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir, a_logdir],
                           print_cmd=False)

        cmd = "cp %s %s" % (os.path.join(a_refdir, self.velmodel), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.stations), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.srffile), a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(a_refdir, self.srcfile), a_indir)
        bband_utils.runprog(cmd)
        for i in range(1, 6):
            cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                             "gp", "s0%d-lf.bbp" % (i)),
                                os.path.join(a_tmpdir,
                                             "%d.s0%d-lf.bbp" % (self.sim_id, i)))
            bband_utils.runprog(cmd)

    def test_bbtoolbox(self):
        """
        Test SDSU BBToolbox code
        """
        bbtool = BBToolbox("", self.velmodel, self.srcfile,
                           self.srffile, self.stations,
                           self.vmodel_name, sim_id=self.sim_id)
        bbtool.run()
        for i in range(1, 6):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "sdsu", "s0%d.bbp" % (i))
            hybfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id),
                                   "%d.s0%d.bbp" % (self.sim_id, i))
            self.failIf(cmp_bbp.cmp_bbp(ref_file, hybfile,
                                        tolerance=0.01) != 0,
                        "output HF BBP %s " % (hybfile) +
                        " file does not match reference hf bbp file %s" %
                        (ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestBBToolbox)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
