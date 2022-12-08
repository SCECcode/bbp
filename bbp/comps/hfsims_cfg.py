#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2022, University of Southern California
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

This config class will encapsulate the configuration parameters
needed to run the GP high-frequency code. Region-specific parameters
can override several of the parameters specified here.
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
        self.HFSIM = "hb_high_v6.1.1"

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
        self.DEFAULT_DX = 1.0
        self.DEFAULT_DY = 1.0

        self.RUPV = -1.0
        self.MEAN_RVFAC = 0.775
        self.RANGE_RVFAC = 0.1
        self.SHAL_RVFAC = 0.6
        self.UNITS = -1
        self.DEFAULT_VSMOHO = 999.9
        self.PATH_DUR_MODEL = 11
        self.DEEP_RVFAC = 0.6
        self.RVSIG = 0.1
        self.C_ZERO = 2.0

        # Extra parameters required by hfsims V5.4
        self.C_ALPHA = -99
        self.FA_SIG1 = 0.0

        # Extra parameter required by hfsims V6.0.3
        # ISPAR_ADJUST = 1 for Western US/Japan
        # ISPAR_ADJUST = 2 for CEUS
        self.ISPAR_ADJUST = 1

        # Extra parameter used by hfsims V6.0.5
        self.SPAR_EXP = 0.5

if __name__ == "__main__":
    print("Test Config Class: %s" % (__file__))
