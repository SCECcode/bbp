#! /usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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
from anderson_gof import AndersonGOF
from station_list import StationList

class TestAndersonGof(unittest.TestCase):
    """
    Unit test for anderson_gof.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.stations = "nr_v13_3_1.stl"
        self.eventname = "NR"
        self.sim_id = int(seqnum.get_seq_num())
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])
        sim_id = self.sim_id

        # Set up paths
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(self.install.A_TMP_DATA_DIR,
                                     str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir, "validations",
                                           "anderson_gof")

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_logdir, a_validation_outdir],
                           print_cmd=False)

        # Copy station and correction files
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "anderson_gof", self.stations),
                            a_indir)
        bband_utils.runprog(cmd, print_cmd=False)

        # Read station list
        slo = StationList(os.path.join(a_indir, self.stations))
        site_list = slo.get_station_list()

        # Loop over stations
        for site in site_list:
            station = site.scode
            src_sims_acc = os.path.join(self.install.A_TEST_REF_DIR,
                                        "anderson_gof", "syn_seis",
                                        "%s.acc.bbp" % (station))
            src_sims_rd50 = os.path.join(self.install.A_TEST_REF_DIR,
                                         "anderson_gof", "syn_seis",
                                         "%s.rd50" % (station))
            dst_sims_acc = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                        (sim_id, station))
            dst_sims_rd50 = os.path.join(a_outdir, "%d.%s.rd50" %
                                         (sim_id, station))
            src_obs_acc = os.path.join(self.install.A_TEST_REF_DIR,
                                       "anderson_gof", "obs_seis",
                                       "%s.bbp" % (station))
            src_obs_rd50 = os.path.join(self.install.A_TEST_REF_DIR,
                                        "anderson_gof", "obs_seis",
                                        "%s.rd50" % (station))
            dst_obs_acc = os.path.join(a_tmpdir_seis, "%s.bbp" %
                                       (station))
            dst_obs_rd50 = os.path.join(a_tmpdir_seis, "%s.rd50" %
                                        (station))

            cmd = "cp %s %s" % (src_sims_acc, dst_sims_acc)
            bband_utils.runprog(cmd, print_cmd=False)

            cmd = "cp %s %s" % (src_sims_rd50, dst_sims_rd50)
            bband_utils.runprog(cmd, print_cmd=False)

            cmd = "cp %s %s" % (src_obs_acc, dst_obs_acc)
            bband_utils.runprog(cmd, print_cmd=False)

            cmd = "cp %s %s" % (src_obs_rd50, dst_obs_rd50)
            bband_utils.runprog(cmd, print_cmd=False)

    def test_anderson_gof(self):
        """
        Run the Anderson GOF test
        """
        gof_obj = AndersonGOF(self.stations, self.eventname, sim_id=self.sim_id)
        gof_obj.run()

        # Check summary GOF file
        ref_sum_file = os.path.join(self.install.A_TEST_REF_DIR, "anderson_gof",
                                    "ref_files",
                                    "gof_anderson.%s.txt" % (self.eventname))
        cal_sum_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                    str(self.sim_id),
                                    "validations",
                                    "anderson_gof",
                                    "%d.gof_anderson.%s.txt" %
                                    (self.sim_id, self.eventname))
        self.assertFalse(cmp_bbp.cmp_files_generic(ref_sum_file, cal_sum_file,
                                                   tolerance=0.005,
                                                   start_col=1) != 0,
                         "GOF Summary file does not match reference file!")

        # Read station list
        slo = StationList(os.path.join(self.install.A_IN_DATA_DIR,
                                       str(self.sim_id), self.stations))
        site_list = slo.get_station_list()

        # Loop over stations
        for site in site_list:
            station = site.scode

            # Check per-station files
            ref_sum_file = os.path.join(self.install.A_TEST_REF_DIR,
                                        "anderson_gof",
                                        "ref_files",
                                        "gof-%s-anderson-%s.txt" %
                                        (self.eventname, station))
            cal_sum_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                        str(self.sim_id),
                                        "validations",
                                        "anderson_gof",
                                        "gof-%s-%d-anderson-%s.txt" %
                                        (self.eventname, self.sim_id,
                                         station))
            self.assertFalse(cmp_bbp.cmp_files_generic(ref_sum_file, cal_sum_file,
                                                       tolerance=0.005,
                                                       start_col=1) != 0,
                             "GOF file for station %s does not match!" %
                             (station))

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestAndersonGof)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
