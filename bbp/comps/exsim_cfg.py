#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
