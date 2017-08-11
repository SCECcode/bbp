#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

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
from station_list import StationList
from install_cfg import InstallCfg
import simplekml
import PlotMap

# Use an extra buffer to plot the region around all stations (in degrees)
BUFFER_LATITUDE = 0.25
BUFFER_LONGITUDE = 0.25

def write_simple_stations(station_file, out_file):
    """
    This function parses the station file and writes a simple
    version with just longitude, latitude, and station code
    """
    stl = StationList(station_file).getStationList()
    fp_out = open(out_file, 'w')
    for stat in stl:
        fp_out.write("%f %f %s\n" % (stat.lon, stat.lat, stat.scode))
    fp_out.flush()
    fp_out.close()

def write_fault_trace(srf_file, out_file):
    """
    This function reads the srf file and outputs a trace file
    """
    points = []
    shallowest = 100.0

    # Figure out SRF file version
    srf_data = open(srf_file, 'r')
    for line in srf_data:
        line = line.strip()
        # Skip blank lines
        if not line:
            continue
        line = int(float(line))
        break

    if line == 1:
        tokens = 8
    elif line == 2:
        tokens = 10
    else:
        bband_utils.ParameterError("Cannot determine SRF file version!")

    for line in srf_data:
        pieces = line.split()
        if len(pieces) == tokens:
            depth = float(pieces[2])
            if depth == shallowest:
                points.append([float(pieces[0]), float(pieces[1])])
            elif depth < shallowest:
                shallowest = depth
                del points[:]
                points.append([float(pieces[0]), float(pieces[1])])
    # Done reading, close file
    srf_data.close()

    # Now, open output file, and write the data
    trace_file = open(out_file, 'w')
    for point in points:
        trace_file.write("%f %f\n" % (point[0], point[1]))
    trace_file.flush()
    trace_file.close()

    # Return trace
    return points

def calculate_fault_edge(lat1, lon1, dist, bearing):
    """
    Given a start point, distance and bearing, calculate the
    destination point
    """
    radius = 6371.0
    to_rad = 0.0174532925

    lat2 = math.asin(math.sin(lat1*to_rad) * math.cos(dist/radius) +
                     math.cos(lat1*to_rad) * math.sin(dist/radius) *
                     math.cos(bearing*to_rad)) / to_rad

    lon2 = (lon1*to_rad +
            math.atan2(math.sin(bearing*to_rad) * math.sin(dist/radius) *
                       math.cos(lat1*to_rad), math.cos(dist/radius) -
                       math.sin(lat1*to_rad) *
                       math.sin(lat2*to_rad))) / to_rad
    return lat2, lon2

def write_simple_trace(a_src_file, out_file):
    """
    This function reads the SRC file and calculates the fault trace
    """
    points = []

    # Read data from SRC file
    cfg_dict = bband_utils.parse_properties(a_src_file)
    if not "fault_length" in cfg_dict:
        raise bband_utils.ParameterError("SRC file missing fault_length!")
    if not "strike" in cfg_dict:
        raise bband_utils.ParameterError("SRC file missing strike!")
    if not "lat_top_center" in cfg_dict:
        raise bband_utils.ParameterError("SRC file missing lat_top_center!")
    if not "lon_top_center" in cfg_dict:
        raise bband_utils.ParameterError("SRC file missing lon_top_center!")
    fault_length = float(cfg_dict["fault_length"])
    strike = float(cfg_dict["strike"])
    lat_top_center = float(cfg_dict["lat_top_center"])
    lon_top_center = float(cfg_dict["lon_top_center"])
    dist = fault_length / 2
    # Calculate 1st edge
    lat1, lon1 = calculate_fault_edge(lat_top_center, lon_top_center,
                                      dist, strike)
    # Reverse direction
    if strike >= 180:
        strike = strike - 180
    else:
        strike = strike + 180
    # Calculate 2nd edge
    lat2, lon2 = calculate_fault_edge(lat_top_center, lon_top_center,
                                      dist, strike)
    points.append([lon1, lat1])
    points.append([lon_top_center, lat_top_center])
    points.append([lon2, lat2])

    # Now, open output file, and write the data
    trace_file = open(out_file, 'w')
    for point in points:
        trace_file.write("%f %f\n" % (point[0], point[1]))
    trace_file.flush()
    trace_file.close()
    # Save trace
    return points

def set_boundaries_from_stations(station_file):
    """
    This function sets the north, south, east, and west boundaries
    of the region we should plot, using the stations' locations in
    the station file
    """
    # Start without anything
    north = None
    south = None
    east = None
    west = None

    # First we read the stations
    stations = StationList(station_file).getStationList()
    # Now go through each one, keeping track of its locations
    for station in stations:
        # If this is the first station, use its location
        if north is None:
            north = station.lat
            south = station.lat
            east = station.lon
            west = station.lon
            # Next station
            continue
        if station.lat > north:
            north = station.lat
        elif station.lat < south:
            south = station.lat
        if station.lon > east:
            east = station.lon
        elif station.lon < west:
            west = station.lon

    # Great, now we just add a buffer on each side
    if north < (90 - BUFFER_LATITUDE):
        north = north + BUFFER_LATITUDE
    else:
        north = 90
    if south > (-90 + BUFFER_LATITUDE):
        south = south - BUFFER_LATITUDE
    else:
        south = -90
    if east < (180 - BUFFER_LONGITUDE):
        east = east + BUFFER_LONGITUDE
    else:
        east = 180
    if west > (-180 + BUFFER_LONGITUDE):
        west = west - BUFFER_LONGITUDE
    else:
        west = -180

    return north, south, east, west

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

    def create_kml_output(self, a_station_file, kml_file,
                          hypo_lat=None, hypo_lon=None):
        """
        Creates a kml output file containing all stations and the fault
        """
        kml = simplekml.Kml()
        stl = StationList(a_station_file).getStationList()
        # Add stations first
        for stat in stl:
            kml.newpoint(name=stat.scode, coords=[(stat.lon, stat.lat)])
        # Now add hypocenter
        if hypo_lat is not None and hypo_lon is not None:
            hyp = kml.newpoint(name="Hypocenter", coords=[(hypo_lon, hypo_lat)])
            hyp.style.iconstyle.color = simplekml.Color.red
        # Add fault trace
        if self.trace is not None and len(self.trace) > 0:
            line_str = kml.newlinestring(name="Fault")
            line_str.altitudemode = simplekml.AltitudeMode.clamptoground
            line_str.style.linestyle.color = simplekml.Color.red
            line_str.style.linestyle.width = 5
            line_str.tessellate = 1
            points = []
            for point in self.trace:
                points.append((point[0], point[1], 0))
            line_str.coords = points
        # Save kml file
        kml.save(kml_file)

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

        a_input_file = os.path.join(a_indir, self.input_file)
        a_station_file = os.path.join(a_indir, self.station_file)

        # Define boundaries to plot using the stations in the station file
        (self.north, self.south,
         self.east, self.west) = set_boundaries_from_stations(a_station_file)

        self.log = os.path.join(install.A_OUT_LOG_DIR, str(self.sim_id),
                                "%d.plot_map.log" % (self.sim_id))
        trace_file = "%s.trace" % (a_input_file)
        simple_station_file = "%s.simple" % (a_station_file)
        if self.input_file.endswith(".srf"):
            self.trace = write_fault_trace(a_input_file, trace_file)
        else:
            self.trace = write_simple_trace(a_input_file, trace_file)
        write_simple_stations(a_station_file, simple_station_file)
        map_prefix = os.path.join(a_outdir, "station_map")
        kml_file = os.path.join(a_outdir, "station_map.kml")
        # Get hypo_lon, hypo_lat from src/srf file
        hypo_lon, hypo_lat = fault_utils.calculate_epicenter(a_input_file)

        # Write the kml file
        self.create_kml_output(a_station_file, kml_file,
                               hypo_lat=hypo_lat, hypo_lon=hypo_lon)

        # Matplotlib
        plottitle = 'Fault Trace with Stations'
        plotregion = [self.west, self.east,
                      self.south, self.north]
        topo = os.path.join(install.A_PLOT_DATA_DIR, 'calTopo18.bf')
        coastal = os.path.join(install.A_PLOT_DATA_DIR, 'gshhs_h.txt')
        border = os.path.join(install.A_PLOT_DATA_DIR, 'wdb_borders_h.txt')
        plotter = PlotMap.PlotMap()
        plotter.plot(plottitle, plotregion, topo,
                     coastal, border, trace_file,
                     simple_station_file, map_prefix,
                     hypo_lat=hypo_lat, hypo_lon=hypo_lon)

        print("Plot MAP Completed".center(80, '-'))

if __name__ == '__main__':
    INPUT_FILE = sys.argv[1]
    STATION_FILE = sys.argv[2]
    PLOT_MAP = Plot_Map(INPUT_FILE, STATION_FILE, int(sys.argv[3]))
    PLOT_MAP.run()
