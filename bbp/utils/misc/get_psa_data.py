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

This module scans rd50 files and compiles information about a specific period
"""

from __future__ import division, print_function

# Import Python modules
import os
import sys

def main():
    """
    Writes PSA information for desired period
    """
    if len(sys.argv) != 4:
        print("Usage: %s station_list input_folder period" %
              (sys.argv[0]))
        sys.exit(-1)

    station_list_input = sys.argv[1]
    input_folder = sys.argv[2]
    period = float(sys.argv[3])

    # Open input file
    input_file = open(station_list_input, 'r')

    # Process each station in the list
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue

        # Got a station
        tokens = line.split()
        if len(tokens) < 3:
            # Skip this line
            continue
        station_lon = tokens[0]
        station_lat = tokens[1]
        station_name = tokens[2]

        # Create PSA input filename
        psa_filename = os.path.join(input_folder, "%s.rd50" % (station_name))

        # Open RotD50 file
        psa_file = open(psa_filename, 'r')

        # Find period
        found_period = False
        for psa_line in psa_file:
            # Skip empty lines and comments
            psa_line.strip()
            if not psa_line:
                continue
            if psa_line.startswith("#") or psa_line.startswith("%"):
                continue

            # Should be a line with data
            psa_tokens = psa_line.split()
            psa_tokens = [psa_token.strip() for psa_token in psa_tokens]
            psa_period = float(psa_tokens[0])

            # Check if this is the period we want
            if psa_period != period:
                continue

            found_period = True
            psa_data = float(psa_tokens[3])
            break

        # Done with this file
        psa_file.close()

        if found_period:
            print("%s %s %.7f" % (station_lat, station_lon, psa_data))

    # Close files
    input_file.close()

if __name__ == "__main__":
    main()
