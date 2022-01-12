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

This module takes care of building a workflow using either user
choices interactively, or an option file containing all needed
parameters.

This code calls bbp2peer for each file in a directory that matches
the prefix and suffix provided. It is used to batch convert observations
files to PEER format to be used in the BBP as recorded data for simulation
comparisons.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import argparse

# Import Broadband modules
import bbp_formatter

def parse_arguments():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Create a BBP simulation "
                                     "from a set of seismograms.")
    parser.add_argument("--input_dir", "-i", dest="input_dir",
                        required=True, help="input directory")
    parser.add_argument("--output_dir", "-o", dest="output_dir",
                        required=True, help="output directory")
    parser.add_argument("--prefix", "-p", dest="prefix",
                        default="",
                        help="prefix for input files")
    parser.add_argument("--suffix", "-s", dest="suffix",
                        default="",
                        help="suffix for input files")
    args = parser.parse_args()

    return args

def read_bbp_header(input_file):
    """
    Reads BBP input header looking for station metadata for writing station list
    """
    station_info = {}

    input_bbp_file = open(input_file, 'r')
    for line in input_bbp_file:
        line = line.strip()
        if line.startswith("#") or line.startswith("%"):
            pieces = line.split()
            if len(pieces) != 3:
                continue
            if pieces[1] == "Station=":
                station_info['station'] = pieces[2]
            elif pieces[1] == "padding=":
                station_info['padding'] = int(pieces[2])
            elif pieces[1] == "hp=":
                station_info['hp'] = float(pieces[2])
            elif pieces[1] == "lp=":
                station_info['lp'] = float(pieces[2])
            elif pieces[1] == "lat=":
                lat = pieces[2]
                if lat[-1].upper() == "N" or lat[-1].upper() == "S":
                    if lat[-1].upper() == "N":
                        sign = 1.0
                    else:
                        sign = -1.0
                    lat = float(lat[:-1]) * sign
                else:
                    lat = float(lat)
                station_info['lat'] = lat
            elif pieces[1] == "lon=":
                lon = pieces[2]
                if lon[-1].upper() == "W" or lon[-1].upper() == "E":
                    if lon[-1].upper() == "E":
                        sign = 1.0
                    else:
                        sign = -1.0
                    lon = float(lon[:-1]) * sign
                else:
                    lon = float(lon)
                station_info['lon'] = lon
        continue
    input_bbp_file.close()

    return station_info

def bbp_2_peer_batch():
    """
    Create a set of PEER files from a directory with acceleration BBP files
    """
    # Get all we need from the command-line
    args = parse_arguments()

    # Get list of matching input files
    files = glob.glob("%s/%s*%s" % (args.input_dir, args.prefix, args.suffix))

    for input_file in sorted(files):
        input_base = os.path.basename(input_file)
        if len(args.suffix):
            suffix_start = input_base.find(args.suffix) - 1
        else:
            suffix_start = None
        station_name = input_base[len(args.prefix):suffix_start]
        # print("Processing: %s ..." % (station_name))
        peer_n = os.path.join(args.output_dir, "%s_N.acc" % (station_name))
        peer_e = os.path.join(args.output_dir, "%s_E.acc" % (station_name))
        peer_z = os.path.join(args.output_dir, "%s_Z.acc" % (station_name))
        bbp_formatter.bbp2peer(input_file, peer_n, peer_e, peer_z)
        station_info = read_bbp_header(input_file)
        if 'station' in station_info:
            print("%.5f %.5f  %s  500.00  %.4f  %.4f" %
                  (station_info['lon'], station_info['lat'],
                   station_info['station'], station_info['hp'] * 1.25,
                   station_info['lp'] * 0.8))

if __name__ == '__main__':
    bbp_2_peer_batch()
