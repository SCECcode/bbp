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

This module reads an event file and output a station list in BBP format.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import argparse

def main():
    """
    Script to extrct a station list from an event file
    """
    parser = argparse.ArgumentParser(description="Extracts a station list "
                                     "from an event file.")
    parser.add_argument("-e", "--event", dest="event_name",
                        required=True, help="event name")
    parser.add_argument("-i", "--input", dest="input_file",
                        required=True, help="input event file")
    parser.add_argument("-o", "--output", dest="output_file",
                        required=True, help="output station file")
    parser.add_argument("-l", "--long-names", dest="long_names",
                        action="store_true", default=False,
                        help="use long station names")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    event_name = args.event_name
    long_names = args.long_names

    # Open files
    event_file = open(input_file, 'r')
    station_file = open(output_file, 'w')

    # Write header
    station_file.write("# BBP Station List for %s\n" % (event_name))
    station_file.write("# Lon    Lat     StationId     Vs30(m/s) ")
    station_file.write("HP_Freq(Hz)   LP_Freq(Hz)\n")

    # Skip event_file header
    _ = event_file.readline()

    # Write stations
    for line in event_file:
        line = line.strip()
        line = line.split()
        long_id = line[0]
        short_id = line[1]
        lat = float(line[2])
        lon = float(line[3])
        vs30 = int(float(line[5]))
        fmax = 1.0 / float(line[6])
        fmin = 1.0 / float(line[7])

        if long_names:
            sta_id = long_id
        else:
            sta_id = short_id

        station_file.write("%7.3f %6.3f  %s %5d %5.4f  %5.4f\n" %
                           (lon, lat, sta_id, vs30, fmin, fmax))

    # Close everything
    event_file.close()
    station_file.close()

if __name__ == "__main__":
    main()
