#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
$Id: hfsims_cfg.py 1744 2016-09-13 22:25:41Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
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

class HfsimsCfg(object):
    """
    Define the configuration parameters for the HFSim program
    """

    def __init__(self, a_srcfile=None):
        """
        Set up some parameters for HFSim
        """

        # Parse src file, if given
        if a_srcfile:
            self.CFGDICT = bband_utils.parse_src_file(a_srcfile)
        else:
            self.CFGDICT = {}

        #
        # Name of executable
        #
        self.HFSIM = "hb_high_v5.4.3"

        #
        # Seismic Parameters
        #
        # As per Rob Graves on 29-November-2012:
        # The SITEAMP parameter is a boolean (0 or 1) flag to tell the code to
        # apply impedance amplification factors related to the prescribed
        # velocity model.  It should be set to 1.
        self.SITEAMP = 1
        # The following parameters are set for western US simulations,
        # override in velocity model configuration file
        self.DEFAULT_SDROP = 50
        self.DEFAULT_FCFAC = 0.0
        self.DEFAULT_QFEXP = 0.6
        self.DEFAULT_C0 = 57
        self.DEFAULT_C1 = 34
        self.RAYSET = [2, 1, 2]

        self.TLEN = 102.4
        # The DT in the high frequency simulation will be set to the
        # value below unless the velocity model used specifies a
        # different high frequency DT value via the HF_DT key.
        self.DT = 0.01

        # As per Rob Graves on 8-February-2013: Actually, the FMAX
        # parameter is not used by the code (but is still required in
        # the input list). The use of KAPPA overrides FMAX in the
        # algorithm.  So, it is kind of a "legacy" parameter.
        # As per Rob Graves on 29-November-2012: The FMAX parameter is
        # over-ridden by KAPPA when KAPPA > 0.  Since all of our runs
        # are using KAPPA > 0, it doesn't matter what we set FMAX to.
        self.FMAX = 10.0
        self.KAPPA = 0.04

        # These are default values for the WUS region
        self.DEFAULT_DX = 2.0
        self.DEFAULT_DY = 2.0

        self.RUPV = -1.0
        self.MEAN_RVFAC = 0.8
        self.RANGE_RVFAC = 0.05
        self.SHAL_RVFAC = 0.7
        self.DEFAULT_EXTRA_FCFAC = 0.0
        self.UNITS = -1
        self.DEFAULT_VSMOHO = 999.9
        self.PATH_DUR_MODEL = 0
        self.DEEP_RVFAC = 1.0
        self.RVSIG = 0.0

        # Extra parameters required by hfsims V5.4
        self.C_ALPHA = 0.1
        self.FA_SIG1 = 0.0

if __name__ == "__main__":
    print("Test Config Class: %s" % (__file__))
