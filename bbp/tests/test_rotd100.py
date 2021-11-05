#!/usr/bin/env python
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
import unittest

# Import Broadband modules
import cmp_bbp
import rotd100
import seqnum
import bband_utils
import install_cfg

class TestRotD100(unittest.TestCase):
    """
    Unit test for the rotd100.py module
    """

    def setUp(self):
        """
        Sets up the environment for the test
        """
        self.install = install_cfg.InstallCfg.getInstance()
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
        bband_utils.mkdirs([self.indir, self.tmpdir, self.outdir, self.logdir],
                           print_cmd=False)

    def test_rotd100(self):
        """
        Test the rotd100 module
        """
        # Load configuration, set sim_id
        sim_id = self.sim_id
        # Reference directory
        ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucb")

        r_station_list = "rotd50.stl"
        a_src_station_list = os.path.join(ref_dir, r_station_list)
        a_dst_station_list = os.path.join(self.indir, r_station_list)

        # Copy station list to indir
        cmd = "cp %s %s" % (a_src_station_list, a_dst_station_list)
        bband_utils.runprog(cmd)

        # Copy velocity bbp files to outdir
        for i in range(1, 6):
            src_bbp = os.path.join(ref_dir, "8000%d.vel.bbp" % (i))
            dst_bbp = os.path.join(self.outdir, "%d.8000%d.vel.bbp" %
                                   (sim_id, i))
            cmd = "cp %s %s" % (src_bbp, dst_bbp)
            bband_utils.runprog(cmd)

        rotd100_obj = rotd100.RotD100(r_station_list, sim_id=sim_id)
        rotd100_obj.run()

        # Check results
        for i in range(1, 6):
            ref_rd100 = os.path.join(ref_dir, "8000%d.rd100" % (i))
            ref_rd100_v = os.path.join(ref_dir, "8000%d.rd100.vertical" % (i))
            new_rd100 = os.path.join(self.outdir, "%d.8000%d.rd100" %
                                     (sim_id, i))
            new_rd100_v = os.path.join(self.outdir, "%d.8000%d.rd100.vertical" %
                                       (sim_id, i))
            self.assertFalse(cmp_bbp.cmp_files_generic(ref_rd100,
                                                       new_rd100) != 0,
                             "Output file %s does not match reference file: %s" %
                             (new_rd100, ref_rd100))
            self.assertFalse(cmp_bbp.cmp_files_generic(ref_rd100_v,
                                                       new_rd100_v) != 0,
                             "Output file %s does not match reference file: %s" %
                             (new_rd100_v, ref_rd100_v))

if __name__ == "__main__":
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestRotD100)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
