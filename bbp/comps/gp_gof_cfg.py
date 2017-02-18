#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This module contains configuration parameters for the GOF generation
$Id: gp_gof_cfg.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import sys

class GPGofCfg(object):
    """
    Define the configuration parameters for gof
    """
    def __init__(self):
        self.COMPS = ["090", "000", "avgh"]
        self.COMPS_PSA5 = ["psa5n", "psa5e", "rotd50"]
        self.MIN_CDST = 0
        self.MAX_CDST = 25

if __name__ == "__main__":
    MSEIS = GPGofCfg()
    print("Created Test Config Class: %s" % (sys.argv[0]))
