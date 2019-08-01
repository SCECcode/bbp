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

Program to create a list of stations for use in the Broadband Platform.
"""
from __future__ import division, print_function

# Import Python modules
import sys

def main():
    """
    Get min and max latitude and longitude values from the user, also
    read the step to be used in generating the station list. This code
    can only be used to generate a maximum of 9999 stations.
    """
    if len(sys.argv) < 6:
        print("Usage: %s lat_min lat_max lon_min lon_max step" % sys.argv[0])
        sys.exit(0)

    lat_min = float(sys.argv[1])
    lat_max = float(sys.argv[2])
    lon_min = float(sys.argv[3])
    lon_max = float(sys.argv[4])
    step = float(sys.argv[5])

    station_number = 0
    cur_lat = lat_min
    cur_lon = lon_min

    while cur_lat <= lat_max:
        while cur_lon <= lon_max:
            station_number = station_number + 1
            print("%2.3f    %2.3f    sta%04d   " % (cur_lon,
                                                    cur_lat,
                                                    station_number))
            cur_lon = cur_lon + step
        cur_lat = cur_lat + step
        cur_lon = lon_min

if __name__ == "__main__":
    main()
