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

This module defines the configuration parameters for the ExSim module
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

# Import Broadband modules
import bband_utils

class ExSimCfg(object):
    """
    Define the configuration parameters for the ExSim program
    """

    def rad(self, deg):
        """
        This function converts from degrees to radians
        """
        # From subroutine rad in converter.for
        return deg * math.atan(1.0)/45.0

    def deg(self, rad):
        """
        This function converts from radians to degrees
        """
        # From subroutine deg in converter.for
        return rad * 45 / math.atan(1.0)

    def prop(self, a, b, c, x):
        """
        From the prop subroutine in converter.for
        """
        d = a * x / math.sqrt(a**2+b**2+c**2)
        e = b * x / math.sqrt(a**2+b**2+c**2)
        f = c * x / math.sqrt(a**2+b**2+c**2)

        return (d, e, f)

    def convert_to_exsim(self):
        """
        This function converts a few BBP SRC parameters into the
        format expected by ExSim
        """
        # From subroutine converter in converter.for
        self.CFGDICT["hypo_along_stk"] = (self.CFGDICT["hypo_along_stk"] +
                                          (self.CFGDICT["fault_length"]/2))
        xlon = self.CFGDICT['lon_top_center']
        xlat = self.CFGDICT['lat_top_center']
        az = self.CFGDICT['strike'] + 180.0
        dis = self.CFGDICT["fault_length"] / 2.0

        # From subroutine azds2cor in converter.for
        alf = dis / 6378
        tet = 90 - xlat
        fi1 = xlon
        teta = self.rad(tet)
        fi = self.rad(fi1)
        azi = self.rad(az)
        p1 = math.cos(alf)
        p2 = math.sin(alf)
        rx = math.sin(teta) * math.cos(fi)
        ry = math.sin(teta) * math.sin(fi)
        rz = math.cos(teta)
        ux = (-math.cos(azi) * math.cos(teta) * math.cos(fi) \
              - math.sin(azi) * math.sin(fi))
        uy = (-math.cos(azi) * math.cos(teta) * math.sin(fi) \
              + math.sin(azi) * math.cos(fi))
        uz = math.cos(azi) * math.sin(teta)
        (rpx, rpy, rpz) = self.prop(rx, ry, rz, p1)
        (upx, upy, upz) = self.prop(ux, uy, uz, p2)

        # From the add subroutine in converter.for
        sx = rpx + upx
        sy = rpy + upy
        sz = rpz + upz

        # From the ctog subroutine in converter.for
        if sy < 0:
            xlon1 = -math.acos(sx / math.sqrt(sx**2 + sy**2))
        else:
            xlon1 = math.acos(sx / math.sqrt(sx**2 + sy**2))
        xlat1 = math.asin(sz / math.sqrt(sx**2 + sy**2 + sz**2))

        self.CFGDICT['lon_top_edge'] = self.deg(xlon1)
        self.CFGDICT['lat_top_edge'] = self.deg(xlat1)

    def calculate_stress(self):
        """
        This function calculates the stress parameter for ExSIM based
        on the depth of the hypocenter
        """
        depth = (self.CFGDICT['depth_to_top'] +
                 self.CFGDICT["hypo_down_dip"] *
                 math.sin(math.radians(self.CFGDICT['dip'])))
        stress = depth * 55 + 175
        return stress

    def __init__(self, a_srcname=None):
        """
        Set up parameters for ExSim
        """
        self.MAX_STATIONS = 300
        self.KAPPA = 0.04
        self.STRESS = 150.0
        self.PARAM_FILE = "EXSIM12.params"
        self.EMPIRICAL_AMPS = "empirical_amps.txt"

        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)
            self.convert_to_exsim()

if __name__ == "__main__":
    EXSIM_CFG = ExSimCfg(sys.argv[1])
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
