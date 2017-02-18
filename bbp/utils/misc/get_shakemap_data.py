#!/usr/bin/env python
"""
Program to calculate the pga and pgv of seismograms from a station
list
$Id: get_shakemap_data.py 1261 2013-10-29 21:53:12Z fsilva $
"""

# Import Python modules
import os
import sys
import glob
import math

def calculate_peak(input_file):
    max_ns = None
    max_ew = None
    my_file = open(input_file, 'r')
    for line in my_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("%") or line.startswith("#"):
            continue
        # This is a good line
        tokens = line.split()
        tokens = [float(token) for token in tokens]
        comp_ns = abs(tokens[1])
        comp_ew = abs(tokens[2])
        if max_ns is None:
            max_ns = comp_ns
            max_ew = comp_ew
            continue
        if comp_ns > max_ns:
            max_ns = comp_ns
        if comp_ew > max_ew:
            max_ew = comp_ew
    my_file.close()

    # Now calculate the geometric mean
    return math.sqrt(max_ns * max_ew)

def main():
    if len(sys.argv) < 3:
        print "Usage: %s station_list output_dir" % (sys.argv[0])
        sys.exit(0)
    station_list = sys.argv[1]
    outdir = sys.argv[2]

    station_file = open(station_list)
    for line in station_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = line.split()
        lon = float(tokens[0])
        lat = float(tokens[1])
        station = tokens[2]
        vel_file = glob.glob("%s/*.%s.vel.bbp" % (outdir, station))[0]
        acc_file = glob.glob("%s/*.%s.acc.bbp" % (outdir, station))[0]
        vel_peak = calculate_peak(vel_file)
        acc_peak = calculate_peak(acc_file)
        print ("%2.3f   %2.3f   %s   %7.5e   %7.5e" %
               (lon, lat, station, vel_peak, acc_peak))

    station_file.close()

if __name__ == "__main__":
    main()
