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
