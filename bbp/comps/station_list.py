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

Class for managing station lists in BB Platform
"""
from __future__ import division, print_function

# Import Python modules
import sys

# Import Broadband modules
from station import Station

# Sets maximum allowed len for station name, code limits are:
# jbsim: 64 characters
# hfsims: 64 characters
# bbtoolbox: 128 characters
# b_green_99v8: 15 characters
# c_simula_v12: 15 characters
# syn1D: 256 characters
MAX_STATION_NAME_LEN = 15

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
