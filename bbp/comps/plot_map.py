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

This Broadband module is used to create the station map file
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

# Import Broadband modules
import bband_utils
import fault_utils
import plot_utils
from station_list import StationList
from install_cfg import InstallCfg
import PlotMap

class Plot_Map(object):

    def __init__(self, input_file, station_file, sim_id=0):
        """
        Initialize class variables
        """
        self.station_file = station_file
        self.input_file = input_file
        self.sim_id = sim_id
        # Box with the region to plot
        self.north = None
        self.south = None
        self.east = None
        self.west = None
        self.trace = None

    def run(self):
        """
        Generates a map showing the fault with stations
        """
        print("Plot MAP".center(80, '-'))

        if (self.input_file is None or self.input_file == "" or
            (not self.input_file.endswith(".srf") and
            not self.input_file.endswith(".src"))):
            # We need a SRC or SRF file to get the fault geometry
            return

        install = InstallCfg.getInstance()

        a_indir = os.path.join(install.A_IN_DATA_DIR, str(self.sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(self.sim_id))

        # Make sure tmpdir exists
        dirs = [a_indir, a_outdir]
        bband_utils.mkdirs(dirs, print_cmd=False)

        a_input_file = os.path.join(a_indir, self.input_file)
        a_station_file = os.path.join(a_indir, self.station_file)

        # Define boundaries to plot using the stations in the station file,
        # and making sure we include the entire fault plane
        (self.north,
         self.south,
         self.east,
         self.west) = plot_utils.set_boundaries_from_stations(a_station_file,
                                                              a_input_file)

        self.log = os.path.join(install.A_OUT_LOG_DIR, str(self.sim_id),
                                "%d.plot_map.log" % (self.sim_id))
        trace_file = "%s.trace" % (a_input_file)
        simple_station_file = "%s.simple" % (a_station_file)
        if self.input_file.endswith(".srf"):
            self.trace = plot_utils.write_fault_trace(a_input_file, trace_file)
        else:
            self.trace = plot_utils.write_simple_trace(a_input_file, trace_file)
        plot_utils.write_simple_stations(a_station_file, simple_station_file)
        map_prefix = os.path.join(a_outdir, "station_map")

        # Get hypo_lon, hypo_lat from src/srf file
        hypo_coord = {}
        hypo_lon, hypo_lat = fault_utils.calculate_epicenter(a_input_file)
        hypo_coord['lat'] = hypo_lat
        hypo_coord['lon'] = hypo_lon

        # Matplotlib
        plottitle = 'Fault Trace with Stations'
        plotregion = [self.west, self.east,
                      self.south, self.north]
        topo = os.path.join(install.A_PLOT_DATA_DIR, 'calTopo18.bf')
        coastal = os.path.join(install.A_PLOT_DATA_DIR, 'gshhs_h.txt')
        border = os.path.join(install.A_PLOT_DATA_DIR, 'wdb_borders_h.txt')
        PlotMap.plot_station_map(plottitle, plotregion, topo,
                                 coastal, border, trace_file,
                                 simple_station_file, map_prefix,
                                 [hypo_coord])

        print("Plot MAP Completed".center(80, '-'))

if __name__ == '__main__':
    INPUT_FILE = sys.argv[1]
    STATION_FILE = sys.argv[2]
    PLOT_MAP = Plot_Map(INPUT_FILE, STATION_FILE, int(sys.argv[3]))
    PLOT_MAP.run()
