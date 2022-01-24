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
            self.assertFalse(cmp_bbp.cmp_bbp(ref_file, hybfile,
                                             tolerance=0.01) != 0,
                             "output HF BBP %s " % (hybfile) +
                             " file does not match reference hf bbp file %s" %
                             (ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestBBToolbox)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
