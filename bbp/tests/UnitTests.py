#! /usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
TS.addTest(unittest.makeSuite(TestPyNGA))
TS.addTest(unittest.makeSuite(TestVm2vm))
TS.addTest(unittest.makeSuite(TestCC))

# Add Graves & Pitarka tests
TS.addTest(unittest.makeSuite(TestGenslip))
TS.addTest(unittest.makeSuite(TestJbsim))
TS.addTest(unittest.makeSuite(TestHfsims))
TS.addTest(unittest.makeSuite(TestWccSiteamp))
TS.addTest(unittest.makeSuite(TestMatch))

# Add Irikura tests
TS.addTest(unittest.makeSuite(TestGenSRF))
TS.addTest(unittest.makeSuite(TestIrikuraHF))

# Add UCSB tests
TS.addTest(unittest.makeSuite(TestUCFaultUtils))
TS.addTest(unittest.makeSuite(TestUCrmg))
TS.addTest(unittest.makeSuite(TestSyn1D))
TS.addTest(unittest.makeSuite(TestUCSite))

# Add SDSU tests
if sys.platform == 'darwin':
    print("*** Mac OS X detected: skipping SDSU BBToolbox unit test.")
else:
    # Don't add on Mac OS X since test will fail due to raytracer issue
    TS.addTest(unittest.makeSuite(TestBBToolbox))

# Add ExSIM tests
TS.addTest(unittest.makeSuite(TestExsim))

# Add Song tests
TS.addTest(unittest.makeSuite(TestRMG))

# Add Post-Processing tests
TS.addTest(unittest.makeSuite(TestRotD50))
TS.addTest(unittest.makeSuite(TestRotD100))
TS.addTest(unittest.makeSuite(TestGPGof))
TS.addTest(unittest.makeSuite(TestSDSUMOGof))
TS.addTest(unittest.makeSuite(TestAndersonGof))
TS.addTest(unittest.makeSuite(TestRZZ2015))
TS.addTest(unittest.makeSuite(TestAS16))

# Done, run the tests
print("==> Running BBP Unit Tests...")
unittest.TextTestRunner(verbosity=2).run(TS)
