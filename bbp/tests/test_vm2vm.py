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

These are unit tests for the vm2vm program
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

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
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())
