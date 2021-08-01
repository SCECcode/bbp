#! /usr/bin/env python
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
"""
from __future__ import division, print_function

# Import Python modules
import os
import shutil
import unittest

# Import Broadband modules
import cmp_bbp
import bband_utils
import seqnum
from install_cfg import InstallCfg
from irikura_hf import IrikuraHF

class TestIrikuraHF(unittest.TestCase):
    """
    Acceptance Test for irikura_hf.py
    """

    def setUp(self):
        """
        Set up and stage in all input files
        """
        self.install = InstallCfg()
        self.velmodel = "nr02-vs500.fk1d"
        self.srcfile = "whittier_v12_11_0_fs.src"
        self.srffile = "whittier_v12_11_0_fs.srf"
        self.stations = "whittier_v19_02_1_short.stl"
        self.stress_drop = "stress_drop.out"
        self.segments_midpoint = "segments.midpoint.txt"
        self.sim_id = int(seqnum.get_seq_num())

        # Set up paths
        refdir = os.path.join(self.install.A_TEST_REF_DIR, "irikura")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir, a_logdir],
                           print_cmd=False)

        shutil.copy2(os.path.join(refdir, self.velmodel),
                     os.path.join(a_indir, self.velmodel))
        shutil.copy2(os.path.join(refdir, self.stations),
                     os.path.join(a_indir, self.stations))
        shutil.copy2(os.path.join(refdir, self.srffile),
                     os.path.join(a_indir, self.srffile))
        shutil.copy2(os.path.join(refdir, self.srcfile),
                     os.path.join(a_indir, self.srcfile))
        shutil.copy2(os.path.join(refdir, self.stress_drop),
                     os.path.join(a_tmpdir, self.stress_drop))
        shutil.copy2(os.path.join(refdir, self.segments_midpoint),
                     os.path.join(a_tmpdir, self.segments_midpoint))

    def test_irikura_hf(self):
        """
        Test Irikura HF code
        """
        irikura_hf_obj = IrikuraHF(self.srcfile, self.srffile,
                                   self.velmodel, self.stations,
                                   "LABasin500", sim_id=self.sim_id)
        irikura_hf_obj.run()
        for i in range(1, 6):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "irikura", "s%02d-hf.bbp" % (i))
            bbpfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id), "%d.s%02d-hf.bbp" %
                                   (self.sim_id, i))
            self.assertFalse(not cmp_bbp.cmp_bbp(bbpfile, ref_file,
                                                 tolerance=0.005) == 0,
                             "output HF BBP file %s " % (bbpfile) +
                             " does not match reference hf bbp file %s" % (ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestIrikuraHF)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
