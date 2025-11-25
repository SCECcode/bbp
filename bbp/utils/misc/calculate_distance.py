#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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
"""

# Import Python modules
from __future__ import division, print_function
import os
import sys
import argparse

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList
import bband_utils

# Import Pynga and its utilities
import pynga.utils as putils

def parse_arguments():
    """
    Parse command-line arguments for the calculate distance command
    """
    parser = argparse.ArgumentParser(description="Calculate station distance(s) from a rupture.")
    parser.add_argument("--station-list",  dest="station_list",
                        help="station list")
    parser.add_argument("--src-file", "--src", dest="src_file",
                        help="fault description in src format")
    parser.add_argument("--lat", dest="latitude", type=float, default=0.0,
                        help="provides the latitude for the station")
    parser.add_argument("--lon", dest="longitude", type=float, default=0.0,
                        help="provides the longitude for the station")
    parser.add_argument("--output", "-o", dest="output_file",
                        help="output file")
    args = parser.parse_args()

    if args.src_file is None:
        print("ERROR: Must specify src file!")
        sys.exit(-1)

    if args.station_list is None:
        if args.latitude == 0.0 or args.longitude == 0.0:
            print("ERROR: Must specify either station list of lat/lon pair!")
            sys.exit(-1)
    else:
        args.latitude = 0.0
        args.lontitude = 0.0

    return args

def calculate_distance():
    """
    Process a station list (or a single station) and compute the
    distance from a rupture specified by a SRC file.
    """
    # Get all we need from the command-line
    args = parse_arguments()

    src_keys = bband_utils.parse_src_file(args.src_file)
    origin = (src_keys['lon_top_center'],
              src_keys['lat_top_center'])
    dims = (src_keys['fault_length'], src_keys['dlen'],
            src_keys['fault_width'], src_keys['dwid'],
            src_keys['depth_to_top'])
    mech = (src_keys['strike'], src_keys['dip'],
            src_keys['rake'])

    if args.station_list is None:
        # Process single-station
        site_geom = [float(args.longitude), float(args.latitude), 0.0]
        (fault_trace1, up_seis_depth,
         low_seis_depth, ave_dip,
         dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
        rjb, rrup, rx = putils.DistanceToSimpleFaultSurface(site_geom,
                                                            fault_trace1,
                                                            up_seis_depth,
                                                            low_seis_depth,
                                                            ave_dip)
        # Print result
        print("%f" % (rrup))
    else:
        # Process station list
        slo = StationList(args.station_list)
        stations = slo.get_station_list()
        for station in stations:
            site_geom = [float(station.lon), float(station.lat), 0.0]
            (fault_trace1, up_seis_depth,
             low_seis_depth, ave_dip,
             dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
            rjb, rrup, rx = putils.DistanceToSimpleFaultSurface(site_geom,
                                                                fault_trace1,
                                                                up_seis_depth,
                                                                low_seis_depth,
                                                                ave_dip)
            print("%s %.4f %.4f %.4f"  % (station.scode, station.lat, station.lon, rrup))

if __name__ == '__main__':
    calculate_distance()
