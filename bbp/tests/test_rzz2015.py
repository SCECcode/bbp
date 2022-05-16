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
from rzz2015 import RZZ2015
from station_list import StationList

class TestRZZ2015(unittest.TestCase):
    """
    Unit test for rzz2015.py
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
                                           "rzz2015")

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_logdir, a_validation_outdir])

        # Copy station list
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "rzz2015", self.stations),
                            a_indir)
        bband_utils.runprog(cmd)

        # Read station list
        slo = StationList(os.path.join(a_indir, self.stations))
        site_list = slo.getStationList()

        # Loop over stations
        for site in site_list:
            station = site.scode
            src_sims_acc = os.path.join(self.install.A_TEST_REF_DIR,
                                        "rzz2015", "syn_seis",
                                        "%s.acc.bbp" % (station))
            dst_sims_acc = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                        (sim_id, station))
            src_obs_acc = os.path.join(self.install.A_TEST_REF_DIR,
                                       "rzz2015", "obs_seis",
                                       "%s.bbp" % (station))
            dst_obs_acc = os.path.join(a_tmpdir_seis, "%s.bbp" %
                                       (station))

            cmd = "cp %s %s" % (src_sims_acc, dst_sims_acc)
            bband_utils.runprog(cmd)

            cmd = "cp %s %s" % (src_obs_acc, dst_obs_acc)
            bband_utils.runprog(cmd)

    def test_rzz2015(self):
        """
        Run the RZZ2015 GOF test
        """
        rzz_obj = RZZ2015(self.stations, self.eventname, sim_id=self.sim_id)
        rzz_obj.run()

        # Check results
        ref_sum_file = os.path.join(self.install.A_TEST_REF_DIR, "rzz2015",
                                    "ref_files",
                                    "rzz2015.%s.txt" % (self.eventname))
        cal_sum_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                    str(self.sim_id),
                                    "validations",
                                    "rzz2015",
                                    "%d.rzz2015.%s.txt" %
                                    (self.sim_id, self.eventname))
        self.assertFalse(cmp_bbp.cmp_files_generic(ref_sum_file, cal_sum_file,
                                                   tolerance=1.0,
                                                   start_col=1,
                                                   sep=",") != 0,
                         "RZZ2015 Summary file does not match reference file!")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestRZZ2015)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
