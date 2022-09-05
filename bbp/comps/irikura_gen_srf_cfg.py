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

class IrikuraGenSrfCfg(object):
    """
    Define the configuration parameters for the Irikura rupture generator
    """

    def __init__(self, a_srcfiles):
        """
        Sets basic class parameters, then parses a_srcnames for more information.
        """
        self.VS = 3.5
        self.DT = 0.025
        self.DENS = 2.7
        self.VEL_RUP_FRAC = 0.8
        self.VEL_RUP_RANGE = 0.0
        self.GENSRF = "gen_srf"
        self.GENSRFSEGMENT = "gen_srf_segment"
        self.SUMSEG = "sum_seg"
        self.CFGDICT = []
        self.num_srcfiles = len(a_srcfiles)

        for a_srcfile in a_srcfiles:
            self.CFGDICT.append(bband_utils.parse_src_file(a_srcfile))

if __name__ == "__main__":
    IRIKURA_GEN_SRF_CFG = IrikuraGenSrfCfg(sys.argv[1:])
    print("Created Test Config Class: %s" % os.path.basename((sys.argv[0])))
