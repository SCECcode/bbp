#!/usr/bin/env python
"""
Copyright 2010-2020 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Define the configuration parameters for the GP rupture generator
"""
from __future__ import division, print_function

# Import Python modules
import os
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

    def __init__(self, a_srcfiles):
        """
        Sets basic class parameters, then parses a_srcname for more information
        """

        # User defined parms
        self.SLIP_SIGMA = 0.75
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

        # Read src files
        self.CFGDICT = []
        self.num_srcfiles = len(a_srcfiles)
        for a_srcfile in a_srcfiles:
            self.CFGDICT.append(bband_utils.parse_src_file(a_srcfile))

if __name__ == "__main__":
    ME = GenslipCfg(sys.argv[1:])
    print("Created Test Config Class: %s" % os.path.basename(sys.argv[0]))
