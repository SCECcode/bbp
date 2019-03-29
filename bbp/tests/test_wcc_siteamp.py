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
        cmd = "mkdir -p %s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = "mkdir -p %s/%d" % (self.install.A_OUT_LOG_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        cmd = ("cp %s/gp/%s %s/%d/." % (self.install.A_TEST_REF_DIR,
                                        self.stations,
                                        self.install.A_IN_DATA_DIR,
                                        self.sim_id))
        bband_utils.runprog(cmd)
        for i in range(1, 6):
            for freq in self.freqs:
                cmd = ("cp %s/gp/s%02d-%s.bbp %s/%d/%d.s%02d-%s.bbp" %
                       (self.install.A_TEST_REF_DIR, i, freq,
                        self.install.A_TMP_DATA_DIR, self.sim_id,
                        self.sim_id, i, freq))
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
                self.failIf(cmp_bbp.cmp_bbp(ref_file, bbpfile) != 0,
                            "output %s BBP file %s does not match reference %s bbp file %s " %
                            (freq, bbpfile, freq, ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestWccSiteamp)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
