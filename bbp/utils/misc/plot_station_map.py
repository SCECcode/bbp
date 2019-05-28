#!/bin/env python
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
