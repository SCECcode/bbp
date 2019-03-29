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
import bband_utils
import seqnum
from install_cfg import InstallCfg
from hfsims_cfg import HfsimsCfg
from hfsims import Hfsims

class TestHfsims(unittest.TestCase):
    """
    Acceptance Test for hfsims.py
    """

    def setUp(self):
        """
        Set up and stage in all input files
        """
        self.install = InstallCfg()
        self.hfsim_cfg = HfsimsCfg()
        self.velmodel = "nr02-vs500.fk1d"
        self.srcfile = "test_wh.src"
        self.srffile = "m5.89-0.20x0.20_s2379646.srf"
        self.stations = "test_stat.txt"
        self.metadata = "metadata.txt"
        self.sim_id = int(seqnum.get_seq_num())

        # Set up paths
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir, a_logdir],
                           print_cmd=False)

        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "gp", self.velmodel), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "gp", self.stations), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "gp", self.srffile), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "gp", self.srcfile), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)

    def test_hfsims(self):
        """
        Test GP HFSims code
        """
        hfs_obj = Hfsims(self.velmodel, self.srcfile, self.srffile, self.stations,
                         "LABasin500", sim_id=self.sim_id)
        hfs_obj.run()
        for i in range(1, 6):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "gp", "s%02d-hf.bbp" % (i))
            bbpfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id), "%d.s%02d-hf.bbp" %
                                   (self.sim_id, i))
            self.failIf(not cmp_bbp.cmp_bbp(bbpfile, ref_file,
                                            tolerance=0.005) == 0,
                        "output HF BBP file %s " % (bbpfile) +
                        " does not match reference hf bbp file %s" % (ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestHfsims)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
