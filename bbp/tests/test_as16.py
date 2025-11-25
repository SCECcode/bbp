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
import math
import unittest

# Import Broadband modules
import cmp_bbp
import seqnum
import bband_utils
from install_cfg import InstallCfg
from as16 import AS16, calculate_as16
from station_list import StationList

def compare_values(val1, val2, tolerance=0.01):
    """
    Check if two values are within a given tolerance,
    return True if yes, or False if no.
    """
    return math.fabs((val1 - val2) / val1) <= tolerance

class TestAS16(unittest.TestCase):
    """
    Unit test for as16.py
    """

    def setUp(self):
        self.install = InstallCfg()
        self.stations = "nr_v19_06_2.stl"
        self.source = "nr_v14_02_1.src"
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
        a_validation_outdir = os.path.join(a_outdir,
                                           "validations",
                                           "stewart_duration_gmpe")

        # Create directories
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_logdir, a_validation_outdir])

        # Copy station list
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "as16", self.stations),
                            a_indir)
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (os.path.join(self.install.A_TEST_REF_DIR,
                                         "as16", self.source),
                            a_indir)
        bband_utils.runprog(cmd)

    def test_as16(self):
        """
        Run the AS16 GMPE test
        """
        as16_obj = AS16(self.stations, self.source,
                        self.eventname, sim_id=self.sim_id)
        as16_obj.run()

        # Check results
        ref_sum_file = os.path.join(self.install.A_TEST_REF_DIR, "as16",
                                    "as16.%s.txt" % (self.eventname))
        cal_sum_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                    str(self.sim_id),
                                    "validations",
                                    "stewart_duration_gmpe",
                                    "%d.as16.%s.txt" %
                                    (self.sim_id, self.eventname))
        self.assertFalse(cmp_bbp.cmp_files_generic(ref_sum_file, cal_sum_file,
                                                   tolerance=0.005,
                                                   start_col=1,
                                                   sep=",") != 0,
                         "AS16 Summary file does not match reference file!")

    def test_as16_testcases(self):
        """
        Run the AS16 test suite
        """
        write_ref_file = False

        # Input file
        ref_filename = os.path.join(self.install.A_TEST_REF_DIR, "as16",
                                    "as16_testcases.txt")
        output_filename = os.path.join(self.install.A_TEST_REF_DIR, "as16",
                                       "as16_python.txt")
        # Process input file
        ref_file = open(ref_filename, 'r')
        if write_ref_file:
            out_file = open(output_filename, 'w')
            out_file.write("Mag,Rrup (km),Focal Mech,Vs30 (m/s),Z1.0 (m),CJ flag,"
                           "Significant Duration 5-75 (s),tau 5-75,phi 5-75,"
                           "Significant Duration 5-95 (s),tau 5-95,phi 5-95,"
                           "Significant Duration 20-80 (s),tau 20-80,phi 20-80\n")
        for line in ref_file:
            line = line.strip()
            if line.startswith("Mag"):
                # Skip header
                continue
            pieces = line.split(",")
            pieces = [float(piece) for piece in pieces]
            pieces[1] = int(pieces[1])
            pieces[2] = int(pieces[2])
            pieces[3] = int(pieces[3])
            pieces[4] = int(pieces[4])
            pieces[5] = int(pieces[5])
            [sd575, sd595, sd2080,
             tau575, tau595, tau2080,
             phi575, phi595, phi2080] = calculate_as16(pieces[0], pieces[1],
                                                       pieces[2], pieces[3],
                                                       pieces[4], pieces[5])
            if write_ref_file:
                out_file.write("%.1f,%d,%d,%d,%d,%d" %
                               (pieces[0], pieces[1], pieces[2], pieces[3],
                                pieces[4], pieces[5]))
                out_file.write(",%7.8f,%7.8f,%7.8f" % (sd575, tau575, phi575))
                out_file.write(",%7.8f,%7.8f,%7.8f" % (sd595, tau595, phi595))
                out_file.write(",%7.8f,%7.8f,%7.8f" % (sd2080, tau2080, phi2080))
                out_file.write("\n")

            # Compare results
            if (compare_values(pieces[6], sd575) and
                compare_values(pieces[7], tau575) and
                compare_values(pieces[8], phi575) and
                compare_values(pieces[9], sd595) and
                compare_values(pieces[10], tau595) and
                compare_values(pieces[11], phi595) and
                compare_values(pieces[12], sd2080) and
                compare_values(pieces[13], tau2080) and
                compare_values(pieces[14], phi2080)):
                # Results are within allowed tolerance, nothing to do
                pass
            else:
                self.assertFalsee(True, "AS16 results do not match reference file!\n" +
                                  "Inputs: %.1f,%d,%d,%d,%d,%d\n" %
                                  (pieces[0], pieces[1], pieces[2],
                                   pieces[3], pieces[4], pieces[5]) +
                                  "Outputs: %7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f\n" %
                                  (sd575, tau575, phi575, sd595,
                                   tau595, phi595, sd2080, tau2080, phi2080) +
                                  "Refs: %7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f,%7.8f\n" %
                                  (pieces[6], pieces[7], pieces[8], pieces[9],
                                   pieces[10], pieces[11], pieces[12], pieces[13],
                                   pieces[14]))

        # Close files
        ref_file.close()
        if write_ref_file:
            out_file.close()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sim_id = int(sys.argv[1])
        print("Will test with sim_id: %d" % (sim_id))
    else:
        sim_id = None
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAS16)
    RETURN_CODE = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(not RETURN_CODE.wasSuccessful())
