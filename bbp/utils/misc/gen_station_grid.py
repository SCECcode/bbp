#!/usr/bin/env python
"""
Program to create a list of stations for use in the Broadband Platform.
$Id: gen_station_grid.py 1260 2013-10-29 21:49:38Z fsilva $
"""

# Import Python modules
import os
import sys

def main():
    """
    Get min and max latitude and longitude values from the user, also
    read the step to be used in generating the station list. This code
    can only be used to generate a maximum of 9999 stations.
    """
    if len(sys.argv) < 6:
        print "Usage: %s lat_min lat_max lon_min lon_max step" % sys.argv[0]
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
            print "%2.3f    %2.3f    sta%04d   10   " % (cur_lon,
                                                         cur_lat,
                                                         station_number)
            cur_lon = cur_lon + step
        cur_lat = cur_lat + step
        cur_lon = lon_min

if __name__ == "__main__":
    main()
