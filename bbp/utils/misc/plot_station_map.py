#!/bin/env python
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

This module plots station map and fault trace
"""
from __future__ import division, print_function

# Import Python Modules
import os
import sys
import argparse

# Import Broadband modules
import PlotMap
import plot_utils
import fault_utils
from install_cfg import InstallCfg

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Creates a map plot with "
                                     " a fault and stations.")
    parser.add_argument("-o", "--output", dest="outfile", required=True,
                        help="output png file")
    parser.add_argument("--station-list", "-sl",
                        dest="station_list", required=True,
                        help="station list with latitude and longitude")
    parser.add_argument("--title", "-t",
                        dest="plot_title",
                        help="title for plot")
    parser.add_argument('src_files', nargs='+')
    args = parser.parse_args()

    return args

def plot_station_map_main():
    """
    Main function for plotting the station map
    """
    # Parse command-line options
    args = parse_arguments()
    # Copy inputs
    output_file = args.outfile
    station_file = args.station_list
    plot_title = args.plot_title
    src_files = args.src_files
    first_src_file = src_files[0]

    # Set paths
    install = InstallCfg.getInstance()
    topo = os.path.join(install.A_PLOT_DATA_DIR, 'calTopo18.bf')
    coastal = os.path.join(install.A_PLOT_DATA_DIR, 'gshhs_h.txt')
    border = os.path.join(install.A_PLOT_DATA_DIR, 'wdb_borders_h.txt')

    # Define boundaries to plot using the stations in the station file,
    # and making sure we include the entire fault plane
    (north, south,
     east, west) = plot_utils.set_boundaries_from_stations(station_file,
                                                           first_src_file)

    trace_file = os.path.join("/tmp",
                              "%s.trace" %
                              (os.path.basename(first_src_file)))
    simple_station_file = os.path.join("/tmp",
                                       "%s.simple" %
                                       (os.path.basename(station_file)))
    trace = plot_utils.write_simple_trace(first_src_file, trace_file)
    plot_utils.write_simple_stations(station_file, simple_station_file)

    # Build a hypocenter list
    hypocenters = []
    for src_file in src_files:
        # Get hypo_lon, hypo_lat from src files
        hypo_lon, hypo_lat = fault_utils.calculate_epicenter(src_file)
        hypo_coord = {}
        hypo_coord['lat'] = hypo_lat
        hypo_coord['lon'] = hypo_lon
        hypocenters.append(hypo_coord)

    # Set plot title
    if plot_title is None:
        plot_title = 'Fault Trace with Stations'
    plot_region = [west, east, south, north]

    # Matplotlib
    PlotMap.plot_station_map(plot_title, plot_region, topo,
                             coastal, border, trace_file,
                             simple_station_file,
                             os.path.splitext(output_file)[0],
                             hypocenters)

    # Delete intermediate files
    os.remove(trace_file)
    os.remove(simple_station_file)

# ============================ MAIN ==============================
if __name__ == "__main__":
    plot_station_map_main()
# end of main program
