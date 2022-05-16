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
from exsim import ExSim
from install_cfg import InstallCfg
import velocity_models as vmodels

class TestExsim(unittest.TestCase):
    """
    Unit test for exsim.py
    """

    def setUp(self):
        """
        Cope needed files to run the test
        """
        self.vmodel_name = "LABasin500"
        self.sim_id = int(seqnum.get_seq_num())
        self.install = InstallCfg()
        self.vmodel_obj = vmodels.get_velocity_model_by_name(self.vmodel_name)

        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir], print_cmd=False)

        # Copy needed files

        # src file
        r_src_file = "nr_v14_02_1.src"
        src_file = os.path.join(self.install.A_TEST_REF_DIR, "uwo", r_src_file)
        self.src_file = os.path.join(indir, r_src_file)
        cmd = "cp %s %s" % (src_file, self.src_file)
        bband_utils.runprog(cmd)

        # exsim param template file
        vmodel_params = self.vmodel_obj.get_codebase_params('exsim')
        self.assertFalse('GENERIC_PARAM' not in vmodel_params)
        r_param_template = vmodel_params['GENERIC_PARAM']

        self.assertFalse(r_param_template == "" or r_param_template is None)
        param_template = os.path.join(self.vmodel_obj.base_dir,
                                      r_param_template)
        # r_param_template is relative to the velocity model basedir,
        # get only basename
        r_param_template = os.path.basename(r_param_template)
        self.param_template = os.path.join(indir, r_param_template)
        cmd = "cp %s %s" % (param_template, self.param_template)
        bband_utils.runprog(cmd)

        # station file
        r_stations = "nr_v19_02_1_one_station.stl"
        stations = os.path.join(self.install.A_TEST_REF_DIR, "uwo", r_stations)
        self.stations = os.path.join(indir, r_stations)
        cmd = "cp %s %s" % (stations, self.stations)
        bband_utils.runprog(cmd)

    def test_exsim(self):
        """
        Run ExSIM module
        """
        exsim_obj = ExSim(self.src_file, self.param_template, self.stations,
                          self.vmodel_name, sim_id=self.sim_id)
        exsim_obj.run()

        for i in range(1, 2):
            ref_file = os.path.join(self.install.A_TEST_REF_DIR, "uwo",
                                    "s%02d.acc.bbp" % (i))
            bbp_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                    str(self.sim_id),
                                    "%d.s%02d.acc.bbp" % (self.sim_id, i))
            self.assertFalse(cmp_bbp.cmp_bbp(ref_file, bbp_file) != 0,
                             "ExSIM output file does not match reference file!")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestExsim)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
