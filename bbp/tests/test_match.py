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
from match_cfg import MatchCfg
from match import Match

class TestMatch(unittest.TestCase):
    """
    Acceptance Test for match.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.match_cfg = MatchCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.stations = "test_stat.txt"
        self.vmodel_name = "LABasin500"
        self.sim_id = int(seqnum.get_seq_num())
        self.freqs = ['lf', 'hf']
        cmd = "mkdir -p %s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_LOG_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                       self.stations,
                                       self.install.A_IN_DATA_DIR,
                                       self.sim_id)
        bband_utils.runprog(cmd)
        for i in range(1, 6):
            for freq in self.freqs:
                cmd = ("cp %s/gp/s%02d-%s-site.bbp %s/%d/%d.s%02d-%s.acc.bbp" %
                       (self.install.A_TEST_REF_DIR, i, freq,
                        self.install.A_TMP_DATA_DIR, self.sim_id,
                        self.sim_id, i, freq))
                bband_utils.runprog(cmd)

    def test_match(self):
        """
        Test GP match code
        """
        match_obj = Match(self.stations, self.vmodel_name, sim_id=self.sim_id)
        match_obj.run()

        for i in range(1, 6):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                    "gp",
                                    "s%02d.merged.bbp" % (i))
            bbpfile = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id),
                                   "%d.s%02d.bbp" % (self.sim_id, i))
            self.failIf(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                        "output merged BBP file %s " % (bbpfile) +
                        " does not match reference merged bbp file %s" %
                        (ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestMatch)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
