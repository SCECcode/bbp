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
import sys
import unittest

# Import Broadband modules
import cmp_bbp
import seqnum
import bband_utils
from install_cfg import InstallCfg
from fas_gof import FASGof

class TestFASGof(unittest.TestCase):
    """
    Unit test for fas_gof.py
    """

    def setUp(self):
        self.install = InstallCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.srcfile = "test_nr.src"
        self.stations = "test_stat_fas.stl"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])

        # Set up paths
        refdir = os.path.join(self.install.A_TEST_REF_DIR, "fas")
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_outdir_fas = os.path.join(a_outdir, "FAS")                                    
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir,
                            a_outdir_fas, a_logdir],
                           print_cmd=False)
        # Copy stations
        cmd = "cp %s %s" % (os.path.join(refdir, self.stations), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)

        # Copy src file
        cmd = "cp %s %s" % (os.path.join(refdir, self.srcfile), a_indir)
        bband_utils.runprog(cmd, print_cmd=False)

        for i in range(11, 14):
            # Copy sample calculated and observed seismograms
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "s%02d.fas.acc.bbp" % (i)),
                                os.path.join(a_outdir,
                                             "%d.s%02d.acc.bbp" %
                                             (self.sim_id, i)))
            bband_utils.runprog(cmd, print_cmd=False)
            # Copy FAS output for both synthetic and observed data
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "s%02d.smc8.smooth.fs.col" % (i)),
                                os.path.join(a_outdir_fas,
                                             "%d.s%02d.smc8.smooth.fs.col" %
                                             (self.sim_id, i)))
            bband_utils.runprog(cmd, print_cmd=False)
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "obs.s%02d.smc8.smooth.fs.col" % (i)),
                                os.path.join(a_outdir_fas,
                                             "obs.s%02d.smc8.smooth.fs.col" % (i)))
            bband_utils.runprog(cmd, print_cmd=False)

    def test_fas_gof(self):
        """
        Test FAS GOF Code
        """
        fas_gof_obj = FASGof(self.srcfile, self.stations,
                             "NR", 'gp', 25, sim_id=self.sim_id)
        fas_gof_obj.run()

        resid_ref_file = os.path.join(self.install.A_TEST_REF_DIR,
                                      "fas", "FASGoF", "NR.fas-resid.txt")
        resid_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                  str(self.sim_id),
                                  "NR-%d.fas-resid.txt" % (self.sim_id))
        self.assertFalse(cmp_bbp.cmp_resid(resid_ref_file,
                                           resid_file,
                                           tolerance=0.005) != 0,
                         "output resid file %s does not match reference resid file %s" %
                         (resid_file, resid_ref_file))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        
        for comp in ['fash1', 'fash2', 'seas']:
            bias_ref_file = os.path.join(self.install.A_TEST_REF_DIR, "fas",
                                         "FASGoF", "NR_r0-25-fas-%s.bias" % (comp))
            m90_ref_file = os.path.join(self.install.A_TEST_REF_DIR, "fas",
                                        "FASGoF", "NR_r0-25-fas-%s.m90" % (comp))
            p90_ref_file = os.path.join(self.install.A_TEST_REF_DIR, "fas",
                                        "FASGoF", "NR_r0-25-fas-%s.p90" % (comp))
            sigma_ref_file = os.path.join(self.install.A_TEST_REF_DIR, "fas",
                                          "FASGoF", "NR_r0-25-fas-%s.sigma" % (comp))
            sigma0_ref_file = os.path.join(self.install.A_TEST_REF_DIR, "fas",
                                           "FASGoF", "NR_r0-25-fas-%s.sigma0" % (comp))
            bias_file = os.path.join(a_outdir, "NR-%d_r0-25-fas-%s.bias" % (self.sim_id, comp))
            m90_file = os.path.join(a_outdir, "NR-%d_r0-25-fas-%s.m90" % (self.sim_id, comp))
            p90_file = os.path.join(a_outdir, "NR-%d_r0-25-fas-%s.p90" % (self.sim_id, comp))
            sigma_file = os.path.join(a_outdir, "NR-%d_r0-25-fas-%s.sigma" % (self.sim_id, comp))
            sigma0_file = os.path.join(a_outdir, "NR-%d_r0-25-fas-%s.sigma0" % (self.sim_id, comp))
            self.assertFalse(cmp_bbp.cmp_bias(bias_ref_file, bias_file) != 0,
                             "output bias file %s does not match reference bias file %s" %
                             (bias_file, bias_ref_file))
            self.assertFalse(cmp_bbp.cmp_bias(m90_ref_file, m90_file) != 0,
                             "output m90 file %s does not match reference m90 file %s" %
                             (m90_file, m90_ref_file))
            self.assertFalse(cmp_bbp.cmp_bias(p90_ref_file, p90_file,
                                              tolerance=0.0025) != 0,
                             "output p90 file %s does not match reference p90 file %s" %
                             (p90_file, p90_ref_file))
            self.assertFalse(cmp_bbp.cmp_bias(sigma_ref_file, sigma_file) != 0,
                             "output sigma file %s does not match reference sigma file %s" %
                             (sigma_file, sigma_ref_file))
            self.assertFalse(cmp_bbp.cmp_bias(sigma0_ref_file, sigma0_file) != 0,
                             "output sigma0 file %s does not match reference sigma0 file %s" %
                             (sigma0_file, sigma0_ref_file))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestFASGof)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
