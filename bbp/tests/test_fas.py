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
from fas import FAS

class TestFAS(unittest.TestCase):
    """
    Unit test for fas.py
    """

    def setUp(self):
        self.install = InstallCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.stations = "test_stat_fas.stl"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])

        # Set up paths
        refdir = os.path.join(self.install.A_TEST_REF_DIR, "fas")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_tmpdir_seis = os.path.join(a_tmpdir, "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_outdir_fas = os.path.join(a_outdir, "FAS")                                    
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_outdir_fas, a_logdir],
                           print_cmd=False)
        # Copy stations
        cmd = "cp %s %s" % (os.path.join(refdir, self.stations), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)

        for i in range(11, 14):
            # Copy sample calculated and observed seismograms
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "s%02d.fas.acc.bbp" % (i)),
                                os.path.join(a_outdir,
                                             "%d.s%02d.acc.bbp" %
                                             (self.sim_id, i)))
            bband_utils.runprog(cmd, print_cmd=False)
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "obs.s%02d.fas.acc.bbp" % (i)),
                                os.path.join(a_tmpdir_seis,
                                             "s%02d.bbp" % (i)))
            bband_utils.runprog(cmd, print_cmd=False)

    def test_fas(self):
        """
        Test FAS Code
        """
        fas_obj = FAS(self.stations, "NR", sim_id=self.sim_id)
        fas_obj.run()

        # Check output
        for i in range(11, 14):
            calc_fas_ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                             "fas", "s%02d.smc8.smooth.fs.col" % (i))
            calc_fas_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                         str(self.sim_id), "FAS",
                                         "%d.s%02d.smc8.smooth.fs.col" % (self.sim_id, i))
            self.assertFalse(cmp_bbp.cmp_fas(calc_fas_file,
                                             calc_fas_ref_file,
                                             tolerance=0.005) != 0,
                             "output fas file %s does not match reference fas file %s" %
                             (calc_fas_file, calc_fas_ref_file))
            
            obs_fas_ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                            "fas", "obs.s%02d.smc8.smooth.fs.col" % (i))
            obs_fas_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                        str(self.sim_id), "FAS",
                                        "obs.s%02d.smc8.smooth.fs.col" % (i))
            self.assertFalse(cmp_bbp.cmp_fas(obs_fas_file,
                                             obs_fas_ref_file,
                                             tolerance=0.005) != 0,
                             "obs fas file %s does not match reference obs fas file %s" %
                             (obs_fas_file, obs_fas_ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestFAS)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
