#!/usr/bin/env python
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

This module compiles results from a single simulation and creates
a csv file containing information about all stations.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import argparse

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList
import bband_utils
import xml_handler

# Initialize global variables
INSTALL = InstallCfg.getInstance()

def parse_arguments():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Create a BBP simulation "
                                     "from a set of seismograms.")
    parser.add_argument("--sim_id", dest="sim_id",
                        type=int, required=True,
                        help="simulation id")
    parser.add_argument("--src_file", "--src", dest="src_file",
                        required=True, help="src file")
    parser.add_argument("--station_list", "--stl", dest="station_list",
                        required=True, help="station list")
    parser.add_argument("--input_dir", "-i", dest="input_dir",
                        required=True, help="input directory")
    parser.add_argument("--prefix", "-p", dest="prefix",
                        default="",
                        help="prefix for input files")
    parser.add_argument("--suffix", "-s", dest="suffix",
                        default="",
                        help="suffix for input files")
    args = parser.parse_args()

    return args

def read_unit_bbp(filename):
    """
    Get the units from the file's header
    Returns either "m" or "cm"
    """
    units = None

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("units=") > 0:
                units = line.split()[2]
                break
        input_file.close()
    except IOError:
        print("[ERROR]: No such file.")
        sys.exit(-1)

    # Make sure we got something
    if units is None:
        print("[ERROR]: Cannot find units in bbp file!")
        sys.exit(-1)

    # Figure out if we have meters or centimeters
    if units == "cm" or units == "cm/s" or units == "cm/s^2":
        return "cm"
    elif units == "m" or units == "m/s" or units == "m/s^2":
        return "m"

    # Invalid units in this file
    print("[ERROR]: Cannot parse units in bbp file!")
    sys.exit(-1)
# end of read_unit_bbp

def copy_seismograms(args, outdata):
    """
    This function copies seismograms into the outdata directory
    """
    slo = StationList(args.station_list)
    site_list = slo.getStationList()

    for site in site_list:
        station_name = site.scode
        print("Importing station: %s ..." % (station_name))
        basename = "%s%s%s" % (args.prefix, station_name, args.suffix)
        base_vel_file = "%s.vel.bbp" % (basename)
        base_acc_file = "%s.acc.bbp" % (basename)
        in_vel_file = os.path.join(args.input_dir, base_vel_file)
        in_acc_file = os.path.join(args.input_dir, base_acc_file)
        out_vel_file = os.path.join(outdata, "%s.%s.vel.bbp" %
                                    (str(args.sim_id), station_name))
        out_acc_file = os.path.join(outdata, "%s.%s.acc.bbp" %
                                    (str(args.sim_id), station_name))

        # Copy vel and acc files
        print("  --> Reading: %s" % (base_vel_file))
        vel_unit = read_unit_bbp(in_vel_file)
        input_file = open(in_vel_file)
        output_file = open(out_vel_file, 'w')
        output_file.write("# Sim NGAH, stat=%s\n" % (station_name))
        output_file.write("#    time(sec)      N-S(cm/s)      E-W(cm/s)      U-D(cm/s)\n")
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("%") or line.startswith("#"):
                continue
            pieces = [float(piece) for piece in line.split()]
            if vel_unit == "m":
                pieces[1] = pieces[1] * 100
                pieces[2] = pieces[2] * 100
                pieces[3] = pieces[3] * 100
            output_file.write("%5.7e\t%5.7e\t%5.7e\t%5.7e\n" %
                              (pieces[0], pieces[1], pieces[2], pieces[3]))
        input_file.close()
        output_file.close()

        print("  --> Reading: %s" % (base_acc_file))
        acc_unit = read_unit_bbp(in_acc_file)
        input_file = open(in_acc_file)
        output_file = open(out_acc_file, 'w')
        output_file.write("#    time(sec)      N-S(cm/s/s)      E-W(cm/s/s)      U-D(cm/s/s)\n")
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("%") or line.startswith("#"):
                continue
            pieces = [float(piece) for piece in line.split()]
            if acc_unit == "m":
                pieces[1] = pieces[1] * 100
                pieces[2] = pieces[2] * 100
                pieces[3] = pieces[3] * 100
            output_file.write("%5.7e\t%5.7e\t%5.7e\t%5.7e\n" %
                              (pieces[0], pieces[1], pieces[2], pieces[3]))
        input_file.close()
        output_file.close()

def create_bbp_sim():
    """
    Create a BBP simulation from a set of seismograms
    """
    # Get all we need from the command-line
    args = parse_arguments()

    # Create simulation directories
    indata = os.path.join(INSTALL.A_IN_DATA_DIR, str(args.sim_id))
    outdata = os.path.join(INSTALL.A_OUT_DATA_DIR, str(args.sim_id))
    tmpdata = os.path.join(INSTALL.A_TMP_DATA_DIR, str(args.sim_id))
    logdir = os.path.join(INSTALL.A_OUT_LOG_DIR, str(args.sim_id))
    bband_utils.mkdirs([indata, outdata, tmpdata, logdir], print_cmd=False)

    # Collect simulation-wide parameters
    copy_seismograms(args, outdata)

if __name__ == '__main__':
    create_bbp_sim()
