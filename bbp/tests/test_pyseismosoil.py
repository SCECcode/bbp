#! /usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2025, University of Southern California
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
from seismo_soil_cfg import SeismoSoilCfg
from seismo_soil import SeismoSoil

class TestPySeismoSoil(unittest.TestCase):
    """
    Acceptance Test for PySeismoSoil site amplification using both Vs30 and Z1.0. 
    """

    def setUp(self):
        self.install = InstallCfg()
        self.vmodel_name = "LABasin500"  
        self.method = "PySeismoSoil"
        print(f"A_INSTALL_ROOT: {self.install.A_INSTALL_ROOT}")
        os.chdir(self.install.A_INSTALL_ROOT)
        self.stations = "test_stat.txt"
        self.velmodel = "test_velmodel.txt"  # Velocity model file
        self.srcfile = "test_src.txt"  # Source file if needed
        self.sim_id = int(seqnum.get_seq_num())

        # Set up directories
        refdir = os.path.join(self.install.A_TEST_REF_DIR, "PySeismoSoil")
        indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        logdir = os.path.join(self.install.A_OUT_LOG_DIR, str(self.sim_id))
        
        # Create all directories
        bband_utils.mkdirs([indir, tmpdir, outdir, logdir], print_cmd=False)

        # Copy required input files
        cmd = "cp %s %s" % (os.path.join(refdir, self.stations), indir)
        bband_utils.runprog(cmd)
        
        cmd = "cp %s %s" % (os.path.join(refdir, self.velmodel), indir)
        bband_utils.runprog(cmd)

        # Copy acceleration BBP files (PySeismoSoil works with acceleration)
        for i in range(1, 6):  # Assuming 5 test stations
            cmd = "cp %s %s" % (os.path.join(refdir,
                                             "s%02d.acc.bbp" % i),
                                os.path.join(outdir,
                                             "%d.s%02d.acc.bbp" %
                                             (self.sim_id, i)))
            bband_utils.runprog(cmd)

    def test_pyseismosoil(self):
        """
        Test PySeismoSoil site amplification module
        """
        # Create PySeismoSoil object
        pyseismo_obj = SeismoSoil(self.srcfile, self.velmodel, self.method,
                                  self.stations, self.vmodel_name, 
                                  sim_id=self.sim_id, debug=False)
        
        # Run the site amplification
        pyseismo_obj.run()

        # Compare results with reference files
        for i in range(1, 6):  # Test stations s01 to s05
            # Test acceleration BBP file
            ref_acc_file = os.path.join(self.install.A_TEST_REF_DIR,
                                        "pyseismosoil",
                                        "s%02d.acc.pyseismo.bbp" % i)
            output_acc_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                           str(self.sim_id),
                                           "%d.s%02d.acc.bbp" %
                                           (self.sim_id, i))
            
            self.assertFalse(cmp_bbp.cmp_bbp(ref_acc_file, output_acc_file) != 0,
                             "PySeismoSoil acceleration output file %s does not match reference file %s" %
                             (output_acc_file, ref_acc_file))
            
            # Test velocity BBP file (if generated)
            ref_vel_file = os.path.join(self.install.A_TEST_REF_DIR,
                                        "pyseismosoil",
                                        "s%02d.vel.pyseismo.bbp" % i)
            output_vel_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                           str(self.sim_id),
                                           "%d.s%02d.vel.bbp" %
                                           (self.sim_id, i))
            
            if os.path.exists(ref_vel_file) and os.path.exists(output_vel_file):
                self.assertFalse(cmp_bbp.cmp_bbp(ref_vel_file, output_vel_file) != 0,
                                 "PySeismoSoil velocity output file %s does not match reference file %s" %
                                 (output_vel_file, ref_vel_file))

    def tearDown(self):
        """
        Clean up test files
        """
        # Optionally remove test directories
        # bband_utils.cleanup(self.sim_id)
        pass

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestPySeismoSoil)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RETURN_CODE.wasSuccessful())