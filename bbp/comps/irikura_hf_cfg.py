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

class IrikuraHFCfg(object):
    """
    Define the configuration parameters for the Irikura Receipe
    Method 2 HF codes
    """
    cfgdict = {}
    vmodel = {}

    def getval(self, attr):
        try:
            val = self.cfgdict[attr]
        except KeyError:
            print("Invalid Source File - Missing attribute: %s" % (attr))
            print("Exiting")
            sys.exit(1)
        return val

    def parse_src(self, a_srcfile):
        """
        This function calls bband_utils' parse property file function
        to get a dictionary of key, value pairs and then looks for the
        parameters needed.
        """
        self.cfgdict = bband_utils.parse_properties(a_srcfile)

        val = self.getval("magnitude")
        self.MAGNITUDE = float(val)

        val = self.getval("fault_length")
        self.LENGTH = float(val)

        val = self.getval("fault_width")
        self.WIDTH = float(val)

        val = self.getval("depth_to_top")
        self.DEPTH_TO_TOP = float(val)

        val = self.getval("strike")
        self.STRIKE = float(val)

        val = self.getval("rake")
        self.RAKE = float(val)

        val = self.getval("dip")
        self.DIP = float(val)

        val = self.getval("lat_top_center")
        self.LAT_TOP_CENTER = float(val)

        val = self.getval("lon_top_center")
        self.LON_TOP_CENTER = float(val)

        val = self.getval("hypo_along_stk")
        self.HYPO_ALONG_STK = float(val)

        val = self.getval("hypo_down_dip")
        self.HYPO_DOWN_DIP = float(val)

        val = self.getval("seed")
        self.SEED = int(val)

    def parse_velmodel(self, a_velmodel):
        """
        This function parses the velocity model file and stores the
        data in the vmodel dictionary.
        """
        # Initialize velocity model structure
        self.vmodel = {'h': [],
                       'vp': [],
                       'vs': [],
                       'rho': [],
                       'qp': [],
                       'qs': []}

        vel_file = open(a_velmodel, 'r')
        for line in vel_file:
            line = line.strip()
            pieces = line.split()
            if len(pieces) == 3:
                self.nlay = float(pieces[0])
                continue
            # Skip lines without the 6 values
            if len(pieces) != 6:
                continue
            pieces = [float(piece) for piece in pieces]
            self.vmodel['h'].append(pieces[0])
            self.vmodel['vp'].append(pieces[1])
            self.vmodel['vs'].append(pieces[2])
            self.vmodel['rho'].append(pieces[3])
            self.vmodel['qp'].append(pieces[4])
            self.vmodel['qs'].append(pieces[5])
        vel_file.close()

    def __init__(self, a_srcname=None, a_velmodel=None):
        """
        Set up parameters for the Irikura recipe
        """
        if a_srcname and a_velmodel:
            self.parse_src(a_srcname)
            self.parse_velmodel(a_velmodel)

        # Filter parameters
        self.filter_order = 3
        self.filter_flo = 1.0e+10
        self.filter_fhi = 0.2
        self.filter_phase = 0
        self.sdropout = "stress_drop.out"
        self.segments_midpoint = "segments.midpoint.txt"

        self.VEL_RUP_FRAC = 0.8
        self.VEL_RUP_RANGE = 0.00
        self.DXX = 1000.0
        self.DYY = 1000.0

if __name__ == "__main__":
    IRIKURA_HF_CFG = IrikuraHFCfg(sys.argv[1])
    print("Created Test Config Class: %s" % (sys.argv[0]))
