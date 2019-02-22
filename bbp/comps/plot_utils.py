#!/usr/bin/env python
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

This Broadband module is used to create the station map file
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

# Import Broadband modules
import bband_utils
from station_list import StationList

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
    line = 0
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

def calculate_fault_edges(a_src_file):
    """
    Calculates the edges of the fault plane
    """
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

    return lat1, lon1, lat_top_center, lon_top_center, lat2, lon2

def write_simple_trace(a_src_file, out_file):
    """
    This function reads the SRC file and calculates the fault trace
    """
    points = []

    (lat1, lon1, lat_top_center,
     lon_top_center, lat2, lon2) = calculate_fault_edges(a_src_file)

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

def set_boundaries_from_stations(station_file, a_src_file):
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

    # Read fault information
    lat1, lon1, _, _, lat2, lon2 = calculate_fault_edges(a_src_file)

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

    # Make sure fault is there Too
    if min(lat1, lat2) < south:
        south = min(lat1, lat2)
    if max(lat1, lat2) > north:
        north = max(lat1, lat2)
    if min(lon1, lon2) < west:
        west = min(lon1, lon2)
    if max(lon1, lon2) > east:
        east = max(lon1, lon2)

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
