#! /usr/bin/python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

These are unit tests for the vm2vm program
"""
from __future__ import division, print_function

# Import Python modules
import os

# Import BBP modules
import vm2vm
import seqnum
import filecmp
import unittest
import bband_utils
from install_cfg import InstallCfg

class TestVm2vm(unittest.TestCase):
    """
    Unit tests for the velocity model conversion codes
    """
    def setUp(self):
        """
        Copy needed files to run the test
        """
        self.install = InstallCfg()
        self.sim_id = int(seqnum.get_seq_num())
        self.a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        self.a_ucsb_refdir = os.path.join(self.install.A_TEST_REF_DIR, "ucsb")
        self.a_sdsu_refdir = os.path.join(self.install.A_TEST_REF_DIR, "sdsu")
        self.a_gp_refdir = os.path.join(self.install.A_TEST_REF_DIR, "gp")

        #
        # Make sure output directories exist
        #
        bband_utils.mkdirs([self.a_tmpdir])

    def test_vm2vm(self):
        """
        input a GP format file and get out a SDSU format file
        """
        gpfile = os.path.join(self.a_sdsu_refdir, "gp_velmodel.v1d")
        sdsufile = os.path.join(self.a_tmpdir, "gp_velmodel.sdsu")
        sdsuref = os.path.join(self.a_sdsu_refdir, "sdsu_velmodel.ref")
        _ = vm2vm.gpvm2sdsuvm(gpfile, sdsufile)
        errmsg = "Conversion of velmodel from GP to SDSU failed"
        self.assertFalse(filecmp.cmp(sdsuref, sdsufile) == False, errmsg)

    def test_vm2vm_nga(self):
        """
        input a GP format file and get out a SDSU format file
        """
        gpfile = os.path.join(self.a_gp_refdir, "nga_rock1.v1d")
        sdsufile = os.path.join(self.a_tmpdir, "sdsu_nga_rock1.v1d")
        sdsuref = os.path.join(self.a_sdsu_refdir, "sdsu_nga_rock1.v1d")
        _ = vm2vm.gpvm2sdsuvm(gpfile, sdsufile)
        errmsg = "Conversion of velmodel from GP to SDSU failed"
        self.assertFalse(filecmp.cmp(sdsuref, sdsufile) == False, errmsg)

    def test_vm2ucsb(self):
        """
        input a GP format file and get out a UCSB format file
        """
        gpfile = os.path.join(self.a_ucsb_refdir, "gp_velocity_model.txt")
        ofile = os.path.join(self.a_tmpdir, "ucsb_velocity_model.txt")
        ucsbref = os.path.join(self.a_ucsb_refdir, "ucsb_velocity_model.txt")
        _ = vm2vm.gpvm2ucsbvm(gpfile, ofile)
        errmsg = "Conversion of velmodel from GP to UCSB ffailed"
        self.assertFalse(filecmp.cmp(ucsbref, ofile) == False, errmsg)

    def test_vm2ucsb_nga(self):
        """
        input a GP format file and get out a UCSB format file
        """
        gpfile = os.path.join(self.a_gp_refdir, "nga_rock1.v1d")
        ofile = os.path.join(self.a_tmpdir, "ucsb_nga_rock1.v1d")
        ucsbref = os.path.join(self.a_ucsb_refdir, "ucsb_nga_rock1.v1d")
        _ = vm2vm.gpvm2ucsbvm(gpfile, ofile)
        errmsg = "Conversion of velmodel from GP to UCSB ffailed"
        self.assertFalse(filecmp.cmp(ucsbref, ofile) == False, errmsg)

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestVm2vm)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
