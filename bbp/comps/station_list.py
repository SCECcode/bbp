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

Class for managing station lists in BB Platform
"""
from __future__ import division, print_function

# Import Python modules
import sys

# Import Broadband modules
from station import Station

# Sets maximum allowed len for station name, code limits are:
# jbsim: 64 characters
# hfsims: 11 characters
# bbtoolbox: 128 characters
# b_green_99v8: 15 characters
# c_simula_v12: 15 characters
# syn1D: 256 characters
MAX_STATION_NAME_LEN = 11

class StationList(object):
    """
    Input Station List file and serve up stations infor as needed
    """
    def __init__(self, a_station_list=None):
        """
        Pass in the absolute file name to the station list and it well be parsed and
        put into a dictionary accessible with iterators.
        """
        if a_station_list is None:
            print("Error reading station list - Null Station List")
            sys.exit(-1)
        self.a_station_filename = a_station_list

        # Start with empty station list
        self.site_list = []

        # Open file
        try:
            station_file = open(self.a_station_filename, "r")
        except IOError:
            print("Error opening station list file : ", a_station_list)
            sys.exit(-1)

        # Read lines one by one
        for line in station_file:
            if line.startswith("#"):
                continue
            sta = line.split()

            if len(sta) >= 3:
                station = Station()
                station.lon = float(sta[0])
                station.lat = float(sta[1])
                station.scode = sta[2]
                if len(station.scode) > MAX_STATION_NAME_LEN:
                    print("Error: station name %s too long!" % (station.scode))
                    print("Maximum limit is %d!" % (MAX_STATION_NAME_LEN))
                    sys.exit(-1)
                if len(sta) >= 4:
                    station.vs30 = int(float(sta[3]))
                if len(sta) >= 6:
                    # We have lf and hf, make sure they are not zero!
                    if float(sta[4]) <= 0:
                        print("warning: station %s has lf<=0, using 1e-15" %
                              (sta[2]))
                        station.low_freq_corner = 1.0e-15
                    else:
                        station.low_freq_corner = float(sta[4])
                    if float(sta[5]) <= 0:
                        print("warning: station %s has hf<=0, using 1e+15" %
                              (sta[2]))
                        station.high_freq_corner = 1.0e+15
                    else:
                        station.high_freq_corner = float(sta[5])
                self.site_list.append(station)
        # Remember to close the file
        try:
            station_file.close
        except IOError:
            print("Error closing station list file :", a_station_list)
        # Error message if we weren't able to read any stations
        if len(self.site_list) == 0:
            print("No stations read from station file :", a_station_list)
            sys.exit(-1)

    @staticmethod
    def build(stat_list, output_file):
        """
        Writes output_file containing a list of our stations
        """
        fp = open(output_file, 'w')
        for stat in stat_list:
            if stat.vs30 is not None:
                fp.write("%f\t%f\t%s\t%d\t%f\t%f\n" %
                         (stat.lon, stat.lat, stat.scode,
                          stat.vs30, stat.low_freq_corner,
                          stat.high_freq_corner))
            else:
                fp.write("%f\t%f\t%f\n" % (stat.lon, stat.lat, stat.scode))
        fp.flush()
        fp.close()

    def getStationList(self):
        """
        Returns our station list
        """
        return self.site_list

if __name__ == "__main__":
    print("Testing Module: %s" % (sys.argv[0]))
    sl = StationList(sys.argv[1])
    ME = sl.getStationList()
    for x in ME:
        print(x.priority)
    sys.exit(0)
