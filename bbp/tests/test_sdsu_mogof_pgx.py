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
import shutil
import unittest

# Import Broadband modules
import cmp_bbp
import seqnum
import bband_utils
from install_cfg import InstallCfg
from sdsu_mogof import SDSUMOGoF

class Test_SDSU_MOGof_PGX(unittest.TestCase):

    def setUp(self):
        """
        Set up unit test
        """
        self.install = InstallCfg()
        os.chdir(self.install.A_INSTALL_ROOT)
        self.r_velocity = ""
        self.gof_weights = dict()
        self.gof_weights["pga"] = 1.0                          #  Weighting on PGA
        self.gof_weights["pgv"] = 1.0                          #  Weighting on PGV
        self.gof_weights["pgd"] = 1.0                          #  Weighting on PGD
        self.gof_weights["psa"] = 1.0                          #  Weighting on PSA
        self.gof_weights["spectral_Fit"] = 0.0                 #  Weighting on Spectral Fit
        self.gof_weights["cumulative_energy_fit"] = 0.0        #  Weighting on Cumulative Energy Fit
        self.gof_weights["inelastic_elastic_fit"] = 0.0        #  Weighting on Inelastic/Elastic Fit (16)
        self.gof_weights["sepctral_acc"] = 0.0                 #  Weighting on Spec Acc (16)
        self.gof_weights["spec_duration"] = 0.0                #  Weighting on Spec Dur (16)
        self.gof_weights["data_energy_release_duration"] = 0.0 #  Weighting on Data Energy Release Duration (5%-75%)
        self.gof_weights["cross_correlation"] = 0.0            #  Weighting on Cross-Correlation
        self.gof_weights["fourier_spectrum"] = 0.0             #  Weighting on Fourier Spectrum

        self.plot_map = True
        self.r_datadir = os.path.join(self.install.A_TEST_REF_DIR, "sdsu")
        self.r_format = "A"
        self.r_mag = 0.0
        self.r_comparison_label = "Landers"
        self.r_cutoff = 10.0
        self.sim_id = int(seqnum.get_seq_num())
        self.a_ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "sdsu")
        self.a_res_dir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        self.in_data_dir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        self.tmp_data_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        self.out_data_dir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        self.out_log_dir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        os.mkdir(self.in_data_dir)
        os.mkdir(self.tmp_data_dir)
        os.mkdir(self.out_data_dir)
        os.mkdir(self.out_log_dir)

        os.chdir(os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id)))

    def tearDown(self):
        os.chdir(self.install.A_TEST_DIR)

    def test_mogof_pgx(self):
        self.r_format = "A"
        self.r_stations = "landers_stats.stl"
        ref_datadir = os.path.join(self.r_datadir, "MOGof")
        cmd = "cp %s/%s %s/%d/." % (ref_datadir, self.r_stations,
                                    self.install.A_IN_DATA_DIR, self.sim_id)
        bband_utils.runprog(cmd)
        ref_inputs_dir = os.path.join(ref_datadir, "ref_inputs", "landers_syn")
        ref_inputs_obs_dir = os.path.join(ref_datadir, "ref_inputs", "landers_obs")
        work_dir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        file_list = os.listdir(ref_inputs_dir)
        file_list = [ref_inputs_dir + "/" + filename for filename in file_list]
        for f in file_list:
            shutil.copy2(f, '%s/%d.%s' % (work_dir, self.sim_id, os.path.basename(f)))

        site_obj = SDSUMOGoF(self.r_stations, self.gof_weights,
                             self.plot_map, ref_inputs_obs_dir,
                             self.r_format, self.r_mag,
                             self.r_comparison_label, self.r_cutoff,
                             sim_id=self.sim_id)
        site_obj.run()

        # Compare individual output files:
        # Assuming we need to compare all metrics
        ref_datadir = os.path.join(ref_datadir, "ref_results", "Landers-UCSB-PGX")
        test_datadir = self.a_res_dir

        self.gof_weights["pgv"] = 1.0                          #  Weighting on PGV
                                  #  Weighting on PGD
                                  #  Weighting on PSA
                         #  Weighting on Spectral Fit
        self.gof_weights["cumulative_energy_fit"] = 0.0        #  Weighting on Cumulative Energy Fit
        self.gof_weights["inelastic_elastic_fit"] = 0.0        #  Weighting on Inelastic/Elastic Fit (16)
        self.gof_weights["sepctral_acc"] = 0.0                 #  Weighting on Spec Acc (16)
        self.gof_weights["spec_duration"] = 0.0                #  Weighting on Spec Dur (16)
        self.gof_weights["data_energy_release_duration"] = 0.0 #  Weighting on Data Energy Release Duration (5%-75%)
        self.gof_weights["cross_correlation"] = 0.0            #  Weighting on Cross-Correlation
        self.gof_weights["fourier_spectrum"] = 0.0             #  Weighting on Fourier Spectrum

        #PGX
        if self.gof_weights["pga"] == 1.0:
            a_ref_file = os.path.join(ref_datadir, "GOF_PGA.list")
            a_newfile = os.path.join(test_datadir, "GOF_PGA.list")
            errmsg = ("Output file %s does not match reference file %s" %
                      (a_newfile, a_ref_file))
            self.assertFalse(cmp_bbp.cmp_gof(a_ref_file, a_newfile,
                                             3, 7, tolerance=0.035) != 0, errmsg)

        if self.gof_weights["pgv"] == 1.0:
            a_ref_file = os.path.join(ref_datadir, "GOF_PGV.list")
            a_newfile = os.path.join(test_datadir, "GOF_PGV.list")
            errmsg = ("Output file %s does not match reference file %s" %
                      (a_newfile, a_ref_file))
            self.assertFalse(cmp_bbp.cmp_gof(a_ref_file, a_newfile,
                                             3, 7, tolerance=0.035) != 0, errmsg)

        if self.gof_weights["pgd"] == 1.0:
            a_ref_file = os.path.join(ref_datadir, "GOF_PGD.list")
            a_newfile = os.path.join(test_datadir, "GOF_PGD.list")
            errmsg = ("Output file %s does not match reference file %s" %
                      (a_newfile, a_ref_file))
            self.assertFalse(cmp_bbp.cmp_gof(a_ref_file, a_newfile,
                                             3, 7, tolerance=0.035) != 0, errmsg)

        if self.gof_weights["psa"] == 1.0:
            a_ref_file = os.path.join(ref_datadir, "GOF_PSA.list")
            a_newfile = os.path.join(test_datadir, "GOF_PSA.list")
            errmsg = ("Output file %s does not match reference file %s" %
                      (a_newfile, a_ref_file))
            self.assertFalse(cmp_bbp.cmp_gof(a_ref_file, a_newfile,
                                             3, 7, tolerance=0.035) != 0, errmsg)

            a_ref_file = os.path.join(ref_datadir, "GOF.list")
            a_newfile = os.path.join(test_datadir, "GOF.list")

            errmsg = ("Output file %s does not match reference file %s" %
                      (a_newfile, a_ref_file))
            self.assertFalse(cmp_bbp.cmp_gof(a_ref_file, a_newfile,
                                             0, 4, tolerance=0.035) != 0, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(Test_SDSU_MOGof_PGX)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
