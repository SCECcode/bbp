#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
$Id: match_cfg.py 1722 2016-08-26 21:50:46Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

class MatchCfg(object):
    """
    Define the configuration parameters for the Match program
    """
    def __init__(self):
        self.HF_FHI = 1.0
        self.HF_FLO = "1.0e+15"
        self.HF_ORD = 4
        self.HF_TSTART = 0.0
        # This should be the same default value used by the hfsims program
        self.NEW_HFDT = 0.01

        self.LF_FHI = 0.0
        self.LF_FLO = 1.0
        self.LF_ORD = 4
        self.LF_TSTART = 0.0

        self.PHASE = 0
        self.MATCH_METHOD = 2

        self.COMPS = ['000', '090', 'ver']

        self.FILTLIST = "filtmatchlist1"

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename((sys.argv[0])))
