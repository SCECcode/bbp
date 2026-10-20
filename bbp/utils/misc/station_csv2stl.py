#!/bin/env python3
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

This module converts a station file from csv format to stl
"""

from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob

def main():
    """
    Converts station file to stl format
    """
    if len(sys.argv) != 3:
        print("Usage: %s input_station_list.csv output_station_list.stl" %
              (sys.argv[0]))
        sys.exit(-1)

    station_list_csv = sys.argv[1]
    station_list_stl = sys.argv[2]

    input_file = open(station_list_csv, 'r')
    output_file = open(station_list_stl, 'w')

    # Walk the station list
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = line.split(",")
        station_id = tokens[0]
        station_lon = float(tokens[1])
        station_lat = float(tokens[2])
        output_file.write("%f\t%f\t%s\n" %
                          (station_lon, station_lat, station_id))
    input_file.close()
    output_file.close()


if __name__ == "__main__":
    main()


