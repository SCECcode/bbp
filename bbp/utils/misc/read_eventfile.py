#!/usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
    station_file.write("LP_Freq(Hz)   HP_Freq(Hz)\n")

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
        hp = 1.0 / float(line[6])
        lp = 1.0 / float(line[7])

        if long_names:
            sta_id = long_id
        else:
            sta_id = short_id

        station_file.write("%7.3f %6.3f  %s %5d %5.4f  %5.4f\n" %
                           (lon, lat, sta_id, vs30, lp, hp))
    
    # Close everything
    event_file.close()
    station_file.close()
        
if __name__ == "__main__":
    main()
