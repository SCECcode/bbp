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
from jbsim_cfg import JbsimCfg
from jbsim import Jbsim

class TestJbsim(unittest.TestCase):
    """
    Acceptance Test for jbsim.py
    """

    def setUp(self):
        """
        Copy needed files to run the test
        """
        self.velmodel = "nr02-vs500.fk1d"
        self.srffile = "m5.89-0.20x0.20_s2379646.srf"
        self.stations = "one_stat.txt"
        self.vmodel_name = "LABasin500"
        self.sim_id = int(seqnum.get_seq_num())

        self.install = InstallCfg()
        self.jbsim_cfg = JbsimCfg(self.vmodel_name)

        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir], print_cmd=False)

        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.velmodel,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.stations,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd, print_cmd=False)
        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.srffile,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd, print_cmd=False)

    def test_jbsim(self):
        """
        Run Jbsim module
        """
        # i_simID,i_r_velmodel,i_r_srffile,i_r_stations,i_r_metadata
        jbsim_obj = Jbsim(self.velmodel, "", self.srffile, self.stations,
                          self.vmodel_name, sim_id=self.sim_id)
        jbsim_obj.run()
        for i in range(1, 2):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "gp",
                                    "s%02d-lf.bbp" % (i))
            bbpfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id),
                                   "%d.s%02d-lf.bbp" % (self.sim_id, i))
            self.failIf(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                        "output LP BBP "
                        "%s file does not match reference lp bbp file %s" %
                        (bbpfile, ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestJbsim)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
