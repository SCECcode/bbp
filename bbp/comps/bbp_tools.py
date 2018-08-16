#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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

    cmd = ("%s/wcc2bbp " % (INSTALL.A_GP_BIN_DIR) +
           "nsfile=%s ewfile=%s udfile=%s " %
           (nsfile, ewfile, udfile) +
           "wcc2bbp=0 < %s >> /dev/null 2>&1" %
           (bbp_in))
    bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)

def join_files(bbp_out, units="cm/s"):
    """
    Joins the 3 1-component bbp files into bbp_out
    """
    nsfile = PREFIX_PROC + COMPS[0]
    ewfile = PREFIX_PROC + COMPS[1]
    udfile = PREFIX_PROC + COMPS[2]

    cmd = ("%s/wcc2bbp " % (INSTALL.A_GP_BIN_DIR) +
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

def plot(bbp_in, png_out):
    """
    Plots the bbp_in file, creating png_out
    """
    plot_seismograms.plot_seis("n/a", bbp_in, "debug", "n/a", png_out)

def comp(bbp1, bbp2, png_out):
    """
    Plots both bbp1 and bbp2, creating png_out
    """
    plot_seismograms.plot_overlay("n/a", bbp1, bbp2, "bbp1", "bbp2", png_out)

def integrate(acc_bbp_in, vel_bbp_out):
    """
    Generates vel_bbp_out by integrating acc_bbp_in
    """
    split_file(acc_bbp_in)
    for component in COMPS:
        filein = PREFIX + component
        fileout = PREFIX_PROC + component
        cmd = ("%s/integ_diff integ=1 filein=%s fileout=%s" %
               (INSTALL.A_GP_BIN_DIR, filein, fileout))
        bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)
    join_files(vel_bbp_out, "cm/s")
    delete_temp_files()

def diff(vel_bbp_in, acc_bbp_out):
    """
    Generates bbp_out by derivating bbp_in
    """
    split_file(vel_bbp_in)
    for component in COMPS:
        filein = PREFIX + component
        fileout = PREFIX_PROC + component
        cmd = ("%s/integ_diff diff=1 filein=%s fileout=%s" %
               (INSTALL.A_GP_BIN_DIR, filein, fileout))
        bband_utils.runprog(cmd, print_cmd=False, abort_on_error=True)
    join_files(acc_bbp_out, "cm/s/s")
    delete_temp_files()

def usage():
    """
    Prints usage and exit
    """
    print("Usage: %s <command> <options>" % os.path.basename(sys.argv[0]))
    print()
    print("Available commands:")
    print(" plot <bbp_in> <png_out>       - plots bbp_in, generating png_out")
    print(" comp <bbp1> <bbp2> <png_out>  - plots bbp1 and bbp2 into png_out")
    print(" integrate <acc_bbp> <vel_bbp> - acc_bbp -> integration -> vel_bbp")
    print(" diff <vel_bbp> <acc_bbp>      - vel_bbp -> diff -> acc_bbp")
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
    plot(sys.argv[2], sys.argv[3])
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
