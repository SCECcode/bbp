#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
$Id: wcc_siteamp_cfg.py 1727 2016-09-01 20:02:22Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

class WccSiteampCfg(object):
    """
    Define the configuration parameters for the Jbrun program
    """
    def __init__(self):
        self.SITEAMP_MODEL3D = "cb2014"
        self.SITEAMP_MODEL = "bssa2014"
        self.FILTLIST = "filtmatchlist1"
        self.GEN_ROCK_VS = 865
        self.VREF_MAX = 1100
        self.FMIN = 0.1
        self.FMIDBOT = 0.2
        self.FLOWCAP = 0.0
        self.COMPS = ["000", "090", "ver"]

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename(sys.argv[0]))
