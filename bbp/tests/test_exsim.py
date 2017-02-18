#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

These are acceptance tests for the exsim.py
$Id: test_exsim.py 1770 2016-10-10 18:13:22Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
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
        self.vmodel_name = "LABasin"
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
        r_src_file = "nr_v12_11_0_fs.src"
        src_file = os.path.join(self.install.A_TEST_REF_DIR, "uwo", r_src_file)
        self.src_file = os.path.join(indir, r_src_file)
        cmd = "cp %s %s" % (src_file, self.src_file)
        bband_utils.runprog(cmd)

        # exsim param template file
        vmodel_params = self.vmodel_obj.get_codebase_params('exsim')
        self.failIf('GENERIC_PARAM' not in vmodel_params)
        r_param_template = vmodel_params['GENERIC_PARAM']

        self.failIf(r_param_template == "" or r_param_template is None)
        param_template = os.path.join(self.vmodel_obj.base_dir,
                                      r_param_template)
        # r_param_template is relative to the velocity model basedir,
        # get only basename
        r_param_template = os.path.basename(r_param_template)
        self.param_template = os.path.join(indir, r_param_template)
        cmd = "cp %s %s" % (param_template, self.param_template)
        bband_utils.runprog(cmd)

        # station file
        r_stations = "nr_v12_11_2.stl"
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
            self.failIf(cmp_bbp.cmp_bbp(ref_file, bbp_file) != 0,
                        "ExSIM output file does not match reference file!")

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestExsim)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
