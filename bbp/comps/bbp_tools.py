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

Broadband Platform Tools
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
from install_cfg import InstallCfg
import bband_utils
import plot_seismograms

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

COMPS = [".000", ".090", ".ver"]
PREFIX = "bbp_tools-temp"
PREFIX_PROC = "bbp_tools-temp.proc"

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

def split_file(bbp_in):
    """
    Splits bbp_in into 3 1-component files
    """
    nsfile = PREFIX + COMPS[0]
    ewfile = PREFIX + COMPS[1]
    udfile = PREFIX + COMPS[2]

    cmd = ("%s " % (os.path.join(INSTALL.A_GP_BIN_DIR, "wcc2bbp")) +
           "nsfile=%s ewfile=%s udfile=%s " %
           (nsfile, ewfile, udfile) +
           "wcc2bbp=0 < %s >> /dev/null 2>&1" %
           (bbp_in))
    bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)

def join_files(bbp_out, units):
    """
    Joins the 3 single-component bbp files into bbp_out
    """
    nsfile = PREFIX_PROC + COMPS[0]
    ewfile = PREFIX_PROC + COMPS[1]
    udfile = PREFIX_PROC + COMPS[2]

    cmd = ("%s " % (os.path.join(INSTALL.A_GP_BIN_DIR, "wcc2bbp")) +
           "nsfile=%s ewfile=%s udfile=%s " %
           (nsfile, ewfile, udfile) +
           "units=%s " % (units) +
           "wcc2bbp=1 > %s 2>> /dev/null" %
           (bbp_out))
    bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)

def delete_temp_files():
    """
    Deletes the 3 single component temporary files
    """
    for component in COMPS:
        filename = PREFIX + component
        os.unlink(filename)
        filename = PREFIX_PROC + component
        os.unlink(filename)

def plot(bbp_in, png_out, station_id, label):
    """
    Plots the bbp_in file, creating png_out
    """
    units_in = read_unit_bbp(bbp_in)
    if units_in == "cm":
        units = "dis"
    elif units_in == "cm/s":
        units = "vel"
    elif units_in == "cm/s/s":
        units = "acc"
    else:
        print("[ERROR]: Cannot parse units in BBP file!")
        sys.exit(-1)

    plot_seismograms.plot_seis(station_id, bbp_in, label, units, png_out)

def comp(bbp1, bbp2, png_out):
    """
    Plots both bbp1 and bbp2, creating png_out
    """
    plot_seismograms.plot_overlay("n/a", bbp1, bbp2, "bbp1", "bbp2", png_out)

def read_unit_bbp(filename):
    """
    Get the units from the file's header
    Returns either "m" or "cm"
    """
    units = None

    try:
        input_file = open(filename, 'r')
        for line in input_file:
            if line.find("time(sec)") > 0:
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

    # Parse and figure what what we got
    units_start = units.find("(")
    units_end = units.find(")")
    if units_start < 0 or units_end < 0:
        print("[ERROR]: Cannot parse units in bbp file!")
        sys.exit(-1)

    units = units[units_start+1:units_end]

    # Check if we got what we needed
    if units == "cm" or units == "cm/s" or units == "cm/s/s":
        return units

    # Invalid units in this file
    print("[ERROR]: Cannot parse units in bbp file!")
    sys.exit(-1)
# end of read_unit_bbp

def integrate(bbp_in, bbp_out):
    """
    Generates bbp_out by integrating bbp_in
    """
    units_in = read_unit_bbp(bbp_in)
    if units_in == "cm":
        print("[ERROR]: Already have a displacement file!")
        sys.exit(-1)
    if units_in == "cm/s/s":
        units_out = "cm/s"
    elif units_in == "cm/s":
        units_out = "cm"
    else:
        print("[ERROR]: Unknown unit in the input file!")
        sys.exit(-1)

    # Split file to get each component separate
    split_file(bbp_in)

    # Process each component
    for component in COMPS:
        filein = PREFIX + component
        fileout = PREFIX_PROC + component
        cmd = ("%s integ=1 filein=%s fileout=%s" %
               (os.path.join(INSTALL.A_GP_BIN_DIR, "integ_diff"),
                filein, fileout))
        bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)

    # Put 3-component BBP file back together
    join_files(bbp_out, units_out)
    delete_temp_files()

def diff(bbp_in, bbp_out):
    """
    Generates bbp_out by derivating bbp_in
    """
    units_in = read_unit_bbp(bbp_in)
    if units_in == "cm/s/s":
        print("[ERROR]: Already have an acceleration file!")
        sys.exit(-1)
    if units_in == "cm":
        units_out = "cm/s"
    elif units_in == "cm/s":
        units_out = "cm/s/s"
    else:
        print("[ERROR]: Unknown unit in the input file!")
        sys.exit(-1)

    # Split file to get each component separate
    split_file(bbp_in)

    # Process each component
    for component in COMPS:
        filein = PREFIX + component
        fileout = PREFIX_PROC + component
        cmd = ("%s diff=1 filein=%s fileout=%s" %
               (os.path.join(INSTALL.A_GP_BIN_DIR, "integ_diff"),
                filein, fileout))
        bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)

    # Put 3-component BBP file back together
    join_files(bbp_out, units_out)
    delete_temp_files()

def usage():
    """
    Prints usage and exit
    """
    print("Usage: %s <command> <options>" % os.path.basename(sys.argv[0]))
    print()
    print("Available commands:")
    print(" plot <bbp_in> <png_out> [station_id run_label] - plots bbp_in")
    print(" comp <bbp1> <bbp2> <png_out>                   - plots bbp1 and bbp2")
    print(" integrate <bbp_in> <bbp_out>                   - integrated bbp_in")
    print(" diff <bbp_in> <bbp_out>                        - differentiates bbp_in")
    print()
    sys.exit(0)

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

INSTALL = InstallCfg.getInstance()

# Check if at least 1 parameter
if len(sys.argv) < 2:
    usage()
# Get command
CMD = sys.argv[1].lower()
if CMD == "plot":
    if len(sys.argv) < 4:
        usage()
    STATION_ID = "n/a"
    LABEL = "n/a"
    if len(sys.argv) > 4:
        STATION_ID = sys.argv[4]
    if len(sys.argv) > 5:
        LABEL = sys.argv[5]
    plot(sys.argv[2], sys.argv[3], STATION_ID, LABEL)
elif CMD == "comp":
    if len(sys.argv) < 5:
        usage()
    comp(sys.argv[2], sys.argv[3], sys.argv[4])
elif CMD == "integrate" or CMD == "int":
    if len(sys.argv) < 4:
        usage()
    integrate(sys.argv[2], sys.argv[3])
elif CMD == "diff":
    if len(sys.argv) < 4:
        usage()
    diff(sys.argv[2], sys.argv[3])
else:
    print("Unknown command: %s" % CMD)
    usage()
