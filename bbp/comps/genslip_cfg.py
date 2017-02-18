#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Define the configuration parameters for the GP rupture generator
$Id: genslip_cfg.py 1727 2016-09-01 20:02:22Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import sys
import random

# Import Broadband modules
import bband_utils

def calculate_rvfac(mean_rvfac, range_rvfac, seed):
    """
    This function calculates a random rvfac value based on the mean
    and range values, plus a seed to generate a random number
    """
    random.seed(seed)
    rvfac = mean_rvfac + range_rvfac * ((random.random() * 2) - 1)
    return rvfac

class GenslipCfg(object):
    """
    Define the configuration parameters for the GP rupture generator
    """

    def __init__(self, a_srcname=None):
        """
        Sets basic class parameters, then parses a_srcname for more information
        """

        # User defined parms
        self.SLIP_SIGMA = 0.85
        # This is now the default inside genslip-3.3, so don't need to use it
        # self.RAND_RAKE_RANGE = 60

        self.RTDEP = 6.5
        self.RTDEP_RANGE = 1.5
        self.MEAN_RVFAC = 0.8
        self.RANGE_RVFAC = 0.05
        self.SHAL_VRUP = 0.6

        # Default RISETIME_COEF set for western US simulations,
        # override in velocity model config file. This parameter used
        # to be set to 1.6, but was modified by RWG in November 2013
        # when the Rupture Generator was updated to version 3.3. The
        # value was reset to 1.6 for Genslip 5.0.1
        self.RISETIME_COEF = 1.6

        # self.EXTRA_RTFAC = 0.0
        self.RISETIME_FAC = 2
        self.RT_SCALEFAC = 1
        self.RT_RAND = 0

        # As in genslip-3.3, we are using 'Mliu' stype, which is the default
        # self.STYPE = "ucsb"

        # Extra parameters in genslip-3.3, updated for genslip-5.0.1
        self.SLIP_WATER_LEVEL = -1
        self.DEEP_RISETIMEDEP = 17.5
        self.DEEP_RISETIMEDEP_RANGE = 2.5
        self.DEEP_RISETIME_FAC = 2.0

        # Read SRC FILE
        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

if __name__ == "__main__":
    ME = GenslipCfg()
    print("Created Test Config Class: %s" % (sys.argv[0]))
