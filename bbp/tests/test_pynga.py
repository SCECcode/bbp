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
"""
from __future__ import division, print_function

# Import Python modules
import os
import math
import unittest

# Import Broadband modules
import pynga
from install_cfg import InstallCfg

class TestPyNGA(unittest.TestCase):
    """
    Unit Tests for PyNGA
    """
    def test_pynga_ngae(self):
        """
        Test PyNGA NGAE
        """
        install = InstallCfg()
        models = ['PZT11', 'A0811E', 'S03SCVS']
        periods = [0.01, 0.04, 0.1, 0.2, 1.0, 2.0]
        
        a_ref_dir = os.path.join(install.A_TEST_REF_DIR, "gmpe")
        a_ref_inputs = os.path.join(a_ref_dir, "NGAE_inputs.csv")
        a_ref_results = os.path.join(a_ref_dir, "NGAE_results.csv")

        #
        # Test each scenario from the input file and compare
        # the result against what we have in the reference file
        #

        input_file = open(a_ref_inputs, 'r')
        results_file = open(a_ref_results, 'r')

        # Skip header line
        line = input_file.readline()
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            tokens = line.split(',')
            tokens = [float(token) for token in tokens]

            # Mw,rjB,rRup,rX,dip,width,zTop,zHyp,rake,vs30,vsInf,z2p5,z1p0
            mag = tokens[0]
            rjb = tokens[1]
            rrup = tokens[2]
            rx = tokens[3]
            dip = tokens[4]
            width = tokens[5]
            ztop = tokens[6]
            zhyp = tokens[7]
            rake = tokens[8]
            vs30 = tokens[9]
            vsinf = tokens[10]
            z2p5 = tokens[11]
            z1p0 = tokens[12]

            if rrup == 0.0:
                rrup = None
            
            if math.isnan(z2p5):
                z2p5 = None
            if math.isnan(z1p0):
                z1p0 = None
            if math.isnan(rake):
                rake = None
                print(line)
                
            for period in periods:
                for model in models:
                    median = pynga.CENA1(model, mag, rjb, rrup, period)
                    calc_val = median
                    ref_line = results_file.readline()
                    ref_line = ref_line.strip()
                    ref_val = float(ref_line)
                    tolerance = 0.01 * max(abs(calc_val), abs(ref_val))  
                    errmsg = ("Calculated and reference values differ!\n" +
                              "Model: %s Period: %f\n" % (model, period) +
                              "Scenario: %s" % (line))
                    self.failIf(abs(calc_val - ref_val) > tolerance, errmsg)
                    # Write reference file
                    #results_file.write("%.9f\n" % (median))
        input_file.close()
        results_file.close()

    def test_pynga_ngaw1(self):
        """
        Test PyNGA NGAW1
        """
        install = InstallCfg()
        models = ['AS', 'BA', 'CB', 'CY']
        periods = [-1, 0.02, 0.2, 1.0, 3.0]
        
        a_ref_dir = os.path.join(install.A_TEST_REF_DIR, "gmpe")
        a_ref_inputs = os.path.join(a_ref_dir, "NGAW1_inputs.csv")
        a_ref_results = os.path.join(a_ref_dir, "NGAW1_results.csv")

        #
        # Test each scenario from the input file and compare
        # the result against what we have in the reference file
        #

        input_file = open(a_ref_inputs, 'r')
        results_file = open(a_ref_results, 'r')

        # Skip header line
        line = input_file.readline()
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            tokens = line.split(',')
            tokens = [float(token) for token in tokens]

            # Mw,rjB,rRup,rX,dip,width,zTop,zHyp,rake,vs30,vsInf,z2p5,z1p0
            mag = tokens[0]
            rjb = tokens[1]
            rrup = tokens[2]
            rx = tokens[3]
            dip = tokens[4]
            width = tokens[5]
            ztop = tokens[6]
            zhyp = tokens[7]
            rake = tokens[8]
            vs30 = tokens[9]
            vsinf = tokens[10]
            z2p5 = tokens[11]
            z1p0 = tokens[12]

            if rrup == 0.0:
                rrup = None
            
            if math.isnan(z2p5):
                z2p5 = None
            if math.isnan(z1p0):
                z1p0 = None
            if math.isnan(rake):
                rake = None
                print(line)
                
            for period in periods:
                for model in models:
                    median, _, _, _ = pynga.NGA08(model, mag, rjb,
                                                  vs30, period,
                                                  rake=rake, dip=dip,
                                                  W=width, Rrup=rrup,
                                                  Rx=rx, Ztor=ztop,
                                                  Zhypo=zhyp, VsFlag=vsinf,
                                                  Z25=z2p5, Z10=z1p0)
                    calc_val = median[0]
                    ref_line = results_file.readline()
                    ref_line = ref_line.strip()
                    ref_val = float(ref_line)
                    tolerance = 0.01 * max(abs(calc_val), abs(ref_val))  
                    errmsg = ("Calculated and reference values differ!\n" +
                              "Model: %s Period: %f\n" % (model, period) +
                              "Scenario: %s" % (line))
                    self.failIf(abs(calc_val - ref_val) > tolerance, errmsg)
                    # Write reference file
                    #results_file.write("%.9f\n" % (median[0]))
        input_file.close()
        results_file.close()
    
    def test_pynga_ngaw2(self):
        """
        Test PyNGA NGAW2
        """
        install = InstallCfg()
        models = ['ASK', 'BSSA', 'CB', 'CY']
        periods = [-1, 0.02, 0.2, 1.0, 3.0]
        
        a_ref_dir = os.path.join(install.A_TEST_REF_DIR, "gmpe")
        a_ref_inputs = os.path.join(a_ref_dir, "NGAW2_inputs.csv")
        a_ref_results = os.path.join(a_ref_dir, "NGAW2_results.csv")

        #
        # Test each scenario from the input file and compare
        # the result against what we have in the reference file
        #

        input_file = open(a_ref_inputs, 'r')
        results_file = open(a_ref_results, 'r')

        # Skip header line
        line = input_file.readline()
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            tokens = line.split(',')
            tokens = [float(token) for token in tokens]

            # Mw,rjB,rRup,rX,dip,width,zTop,zHyp,rake,vs30,vsInf,z2p5,z1p0
            mag = tokens[0]
            rjb = tokens[1]
            rrup = tokens[2]
            rx = tokens[3]
            dip = tokens[4]
            width = tokens[5]
            ztop = tokens[6]
            zhyp = tokens[7]
            rake = tokens[8]
            vs30 = tokens[9]
            vsinf = tokens[10]
            z2p5 = tokens[11]
            z1p0 = tokens[12]

            if rrup == 0.0:
                rrup = None
            
            if math.isnan(z2p5):
                z2p5 = None
            if math.isnan(z1p0):
                z1p0 = None
            
            for period in periods:
                for model in models:
                    median, _, _, _ = pynga.NGA14(model, mag, rjb,
                                                  vs30, period,
                                                  rake=rake, dip=dip,
                                                  W=width, Rrup=rrup,
                                                  Rx=rx, Ztor=ztop,
                                                  Zhypo=zhyp, VsFlag=vsinf,
                                                  Z25=z2p5, Z10=z1p0)
                    calc_val = median[0]
                    ref_line = results_file.readline()
                    ref_line = ref_line.strip()
                    ref_val = float(ref_line)
                    tolerance = 0.01 * max(abs(calc_val), abs(ref_val))  
                    errmsg = ("Calculated and reference values differ!\n" +
                              "Model: %s Period: %f\n" % (model, period) +
                              "Scenario: %s" % (line))
                    self.failIf(abs(calc_val - ref_val) > tolerance, errmsg)
                    # Write reference file
                    #results_file.write("%.9f\n" % (median[0]))
        input_file.close()
        results_file.close()

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestPyNGA)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
