#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Define the configuration parameters for the GP rupture generator
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import random

# Import Broadband modules
import bband_utils

def calculate_rvfac(mean_rvfac, range_rvfac, seed, count=0):
    """
    This function calculates a random rvfac value based on the mean
    and range values, plus a seed to generate a random number
    """
    random.seed(seed)
    while count > 0:
        count = count - 1
        random.seed(int(random.random() * 100000000))
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
        self.MEAN_RVFAC = 0.775
        self.RANGE_RVFAC = 0.1
        self.SHAL_VRUP = 0.6

        # Default RANGE_FWIDTH_FRAC value (randomization disabled)
        self.RANGE_FWIDTH_FRAC = 0.0

        # Default RISETIME_COEF set for western US simulations,
        # override in velocity model config file. This parameter used
        # to be set to 1.6, but was modified by RWG in November 2013
        # when the Rupture Generator was updated to version 3.3. The
        # value was reset to 1.6 for Genslip 5.0.1, and then set to
        # 2.3 for genslip-5.5.2
        self.RISETIME_COEF = 2.3

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
