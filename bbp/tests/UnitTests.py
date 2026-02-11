#! /usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2026, University of Southern California
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

Top level test suites for BB Platform
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import unittest

# Import Broadband modules
from install_cfg import InstallCfg
from test_bband_core import CoreTestSuite
from test_pynga import TestPyNGA
from test_genslip import TestGenslip
from test_jbsim import TestJbsim
from test_hfsims import TestHfsims
from test_wcc_siteamp import TestWccSiteamp
from test_match import TestMatch
from test_gensrf import TestGenSRF
from test_irikura_hf import TestIrikuraHF
from test_uc_fault_utils import TestUCFaultUtils
from test_vm2vm import TestVm2vm
from test_cc import TestCC
from test_ucrmg import TestUCrmg
from test_syn1d import TestSyn1D
from test_uc_site import TestUCSite
from test_bbtoolbox import TestBBToolbox
from test_exsim import TestExsim
from test_rmg import TestRMG
from test_rotd50 import TestRotD50
from test_rotd100 import TestRotD100
from test_gp_gof import TestGPGof
from test_sdsu_mogof import TestSDSUMOGof
from test_anderson_gof import TestAndersonGof
from test_rzz2015 import TestRZZ2015
from test_as16 import TestAS16
from test_fas import TestFAS
from test_fas_gof import TestFASGof

class Logger(object):
    def __init__(self, filename):
        self.filename = filename
        self.out_fp = open(self.filename, 'w')

    def write(self, string):
        self.out_fp.write(string)

    def flush(self):
        self.out_fp.flush()

    def close(self):
        self.out_fp.flush()
        self.out_fp.close()

# Initialize
INSTALL = InstallCfg.getInstance()
sys.stdout = Logger(os.path.join(INSTALL.A_OUT_LOG_DIR, "unit_tests.log"))
TS = unittest.TestSuite()

# Add broadband platform generic tests
TS.addTests(CoreTestSuite())
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestPyNGA))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestVm2vm))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestCC))

# Add Graves & Pitarka tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestGenslip))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestJbsim))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestHfsims))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestWccSiteamp))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestMatch))

# Add Irikura tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestGenSRF))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestIrikuraHF))

# Add UCSB tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestUCFaultUtils))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestUCrmg))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestSyn1D))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestUCSite))

# Add SDSU tests
if sys.platform == 'darwin':
    print("*** Mac OS X detected: skipping SDSU BBToolbox unit test.")
else:
    # Don't add on Mac OS X since test will fail due to raytracer issue
    TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestBBToolbox))

# Add ExSIM tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestExsim))

# Add Song tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestRMG))

# Add Post-Processing tests
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestRotD50))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestRotD100))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestGPGof))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestFAS))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestFASGof))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestSDSUMOGof))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestAndersonGof))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestRZZ2015))
TS.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestAS16))

# Done, run the tests
print("==> Running BBP Unit Tests...")
RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(TS)
sys.exit(not RETURN_CODE.wasSuccessful())
