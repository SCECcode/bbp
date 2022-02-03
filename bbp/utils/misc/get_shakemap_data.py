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

Program to calculate the pga and pgv of seismograms from a station
list
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import math

def calculate_peak(input_file):
    max_ns = None
    max_ew = None
    my_file = open(input_file, 'r')
    for line in my_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("%") or line.startswith("#"):
            continue
        # This is a good line
        tokens = line.split()
        tokens = [float(token) for token in tokens]
        comp_ns = abs(tokens[1])
        comp_ew = abs(tokens[2])
        if max_ns is None:
            max_ns = comp_ns
            max_ew = comp_ew
            continue
        if comp_ns > max_ns:
            max_ns = comp_ns
        if comp_ew > max_ew:
            max_ew = comp_ew
    my_file.close()

    # Now calculate the geometric mean
    return math.sqrt(max_ns * max_ew)

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: %s station_list output_dir [pga|pgv|pga-log|pga-g|pga-g-log|pgv-log]" %
              (sys.argv[0]))
        sys.exit(-1)
    station_list = sys.argv[1]
    outdir = sys.argv[2]
    output_mode = "all"
    if len(sys.argv) == 4:
        output_mode = sys.argv[3].lower()
    if output_mode not in ["all", "pga", "pgv", "pga-log",
                           "pga-g", "pga-g-log", "pgv-log"]:
        print("Requested output mode is invalid!")
        sys.exit(-1)

    station_file = open(station_list)
    for line in station_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = line.split()
        lon = float(tokens[0])
        lat = float(tokens[1])
        station = tokens[2]

        # Calculate peak values
        vel_file = glob.glob("%s/*.%s.vel.bbp" % (outdir, station))[0]
        acc_file = glob.glob("%s/*.%s.acc.bbp" % (outdir, station))[0]
        vel_peak = calculate_peak(vel_file)
        acc_peak = calculate_peak(acc_file)

        # Convert units if needed
        if output_mode == "pga-log":
            acc_peak = math.log10(acc_peak)
        elif output_mode == "pga-g":
            acc_peak = acc_peak / 980.665
        elif output_mode == "pga-g-log":
            acc_peak = math.log10(acc_peak / 980.665)
        elif output_mode == "pgv-log":
            vel_peak = math.log10(vel_peak)

        # Print output
        if output_mode == "all":
            print("%2.3f   %2.3f   %s   %7.5e   %7.5e" %
                  (lon, lat, station, vel_peak, acc_peak))
        elif output_mode in ["pga", "pga-log", "pga-g", "pga-g-log"]:
            print("%2.3f   %2.3f   %7.5e" % (lon, lat, acc_peak))
        elif output_mode in ["pgv", "pgv-log"]:
            print("%2.3f   %2.3f   %7.5e" % (lon, lat, vel_peak))
        else:
            print("Invalid output mode!")
            sys.exit(-1)

    station_file.close()

if __name__ == "__main__":
    main()
