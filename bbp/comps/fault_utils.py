#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This file includes several utility funcions used by various BBP
modules to get information about the rupture.
"""
from __future__ import division, print_function

# Import Python modules
import os
import math
import uuid
import shutil
import tempfile

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg

def get_magnitude(velfile, srffile, suffix="tmp"):
    """
    Scans the srffile and returns the magnitude of the event
    """
    magfile = os.path.join(tempfile.gettempdir(),
                           "%s_%s" %
                           (str(uuid.uuid4()), suffix))
    install = InstallCfg.getInstance()
    cmd = ("%s velfile=%s < %s 2> %s" %
           (os.path.join(install.A_GP_BIN_DIR, "srf2moment"),
            velfile, srffile, magfile))
    bband_utils.runprog(cmd, False)
    srf2moment_fp = open(magfile, 'r')
    srf2moment_data = srf2moment_fp.readlines()
    srf2moment_fp.close()
    #magnitude on last line
    mag_line = srf2moment_data[len(srf2moment_data) - 4]
    pieces = mag_line.split()
    magnitude = float(pieces[5].split(")")[0])
    cmd = "rm %s" % (magfile)
    bband_utils.runprog(cmd, False)
    return magnitude

def get_hypocenter(srffile, suffix="tmp"):
    """
    Looks up the hypocenter of an event in a srffile
    """
    hypfile = os.path.join(tempfile.gettempdir(),
                           "%s_%s" %
                           (str(uuid.uuid4()), suffix))
    install = InstallCfg.getInstance()
    cmd = ("%s < %s > %s" %
           (os.path.join(install.A_GP_BIN_DIR, "srf_gethypo"),
            srffile, hypfile))
    bband_utils.runprog(cmd)
    srf_hypo_fp = open(hypfile, 'r')
    srf_hypo_data = srf_hypo_fp.readline()
    srf_hypo_fp.close()
    srf_hypo = srf_hypo_data.split()
    hypo = []
    for i in range(0, 3):
        hypo.append(float(srf_hypo[i]))
    cmd = "rm %s" % (hypfile)
    bband_utils.runprog(cmd)
    return hypo

def calculate_hypo_depth(srcfile):
    """
    Calculates the hypocenter depth using the SRC file parameters
    """
    cfgdict = bband_utils.parse_properties(srcfile)

    # Look for the needed keys in the SRC file
    try:
        depth_to_top = cfgdict["depth_to_top"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing DEPTH_TO_TOP parameter!")
    depth_to_top = float(depth_to_top)

    try:
        dip = cfgdict["dip"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing DIP parameter!")
    dip = float(dip)

    try:
        hypo_down_dip = cfgdict["hypo_down_dip"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing HYPO_DOWN_DIP parameter!")
    hypo_down_dip = float(hypo_down_dip)

    # Now, calculate the hypocenter depth
    hypo_depth = depth_to_top + hypo_down_dip * math.sin(math.radians(dip))

    # Done
    return hypo_depth

def calculate_epicenter(input_file):
    """
    This function returns the epicenter of an event using either a SRC
    file or a SRF file to look for the hypocenter location. It uses
    Rob Graves' xy2ll utility to convert the coordinates to lat/lon.
    """
    # If we have a SRF file, we already have a function that does this
    if input_file.endswith(".srf"):
        # Get information from srf file
        hypo_lon, hypo_lat, _ = get_hypocenter(input_file)
        return hypo_lon, hypo_lat

    # If we don't have a SRC file, we should print an error here
    if not input_file.endswith(".src"):
        bband_utils.ParameterError("input file should be a SRC or SRF file!")

    # Ok, we have a SRC file
    # Get information from SRC file
    cfgdict = bband_utils.parse_properties(input_file)

    try:
        strike = cfgdict["strike"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing STRIKE parameter!")
    strike = float(strike)

    try:
        dip = cfgdict["dip"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing DIP parameter!")
    dip = float(dip)

    try:
        hypo_down_dip = cfgdict["hypo_down_dip"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing "
                                   "HYPO_DOWN_DIP parameter!")
    hypo_down_dip = float(hypo_down_dip)

    try:
        hypo_along_strike = cfgdict["hypo_along_stk"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing "
                                   "HYPO_ALONG_STK parameter!")
    hypo_along_strike = float(hypo_along_strike)

    try:
        lat_top_center = cfgdict["lat_top_center"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing "
                                   "LAT_TOP_CENTER parameter!")
    lat_top_center = float(lat_top_center)

    try:
        lon_top_center = cfgdict["lon_top_center"]
    except KeyError:
        bband_utils.ParameterError("SRC file missing "
                                   "LON_TOP_CENTER parameter!")
    lon_top_center = float(lon_top_center)

    # Ok, we have all the parameters that we need!
    hypo_perpendicular_strike = hypo_down_dip * math.cos(math.radians(dip))

    # Now call xy2ll program to convert it to lat/long
    # Create temp directory to avoid any race conditions
    tmpdir = tempfile.mkdtemp(prefix="bbp-")
    hypfile = os.path.join(tmpdir, "src_hypo.tmp")
    install = InstallCfg.getInstance()
    cmd = ('echo "%f %f" | %s mlat=%f mlon=%f xazim=%f > %s' %
           (hypo_along_strike, hypo_perpendicular_strike,
            os.path.join(install.A_GP_BIN_DIR, "xy2ll"),
            lat_top_center, lon_top_center, strike, hypfile))
    bband_utils.runprog(cmd, print_cmd=False)
    src_hypo_fp = open(hypfile, 'r')
    src_hypo_data = src_hypo_fp.readline()
    src_hypo_fp.close()
    src_hypo = [float(val) for val in src_hypo_data.split()]
    # Delete temp directory
    shutil.rmtree(tmpdir)

    # Return calculated lon/lat
    return src_hypo[0], src_hypo[1]
