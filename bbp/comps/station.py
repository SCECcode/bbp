#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Class for managing station info in BB Platform
$Id: station.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import sys

class Station(object):
    """
    Class for storing station information
    """

    lon = None
    lat = None
    scode = None
    vs30 = None
    low_freq_corner = 1.0e-15
    high_freq_corner = 1.0e+15

    # SDSU fields
    x = None
    y = None
    vp = None
    vs = None
    rho = None
    kappa = None

if __name__ == "__main__":
    print("Testing Module: %s" % (sys.argv[0]))
    ME = Station()
    ME.lon = 34.00
    ME.lat = -118.00
    ME.scode = "s003"
    ME.vs30 = 865
    sys.exit(0)
