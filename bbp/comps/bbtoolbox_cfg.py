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

This module defines the configuration parameters for the BBToolbox script
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import cc
import bband_utils

class BBToolboxCfg(object):
    """
    Define the configuration parameters for the SDSU BBToolbox program
    """
    cfgdict = {}

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
        This function calls bband_utils's parse property file function
        to get a dictionary of key, value pairs and then looks for a
        the parameters needed by bbtoolbox
        """
        self.cfgdict = bband_utils.parse_properties(a_srcfile)

        val = self.getval("depth_to_top")
        self.DEPTH_TO_TOP = float(val)

        val = self.getval("fault_length")
        self.LENGTH = float(val)

        val = self.getval("dip")
        self.DIP = float(val)

        val = self.getval("rake")
        self.RAKE = float(val)

        val = self.getval("hypo_along_stk")
        self.HYPO_ALONG_STK = float(val)

        val = self.getval("hypo_down_dip")
        self.HYPO_DOWN_DIP = float(val)

        val = self.getval("magnitude")
        self.MAG = float(val)

        val = self.getval("seed")
        self.SEED = int(float(val))

        # Now look for the optional grid parameters
        if 'grid_x' in self.cfgdict:
            self.grid_x = float(self.getval("grid_x"))
        if 'grid_y' in self.cfgdict:
            self.grid_y = float(self.getval("grid_y"))
        if 'grid_z' in self.cfgdict:
            self.grid_z = float(self.getval("grid_z"))

        #
        # Read parameters out of the source file to obtain parameters
        # needed by the BBcoda codes
        #
        fcodes = cc.find_fx_fy_fz(self.HYPO_ALONG_STK,
                                  self.LENGTH,
                                  self.DIP,
                                  self.HYPO_DOWN_DIP,
                                  self.DEPTH_TO_TOP)
        self.fsx = fcodes[0]
        self.fsy = fcodes[1]
        self.fsz = fcodes[2]
        #print ("ETH conversion from hypalongstk: "
        #       "%f flength: %f dip: %f hypdowndip: %f depthtotop: %f\n" %
        #       (self.HYPO_ALONG_STK,
        #        self.LENGTH,
        #        self.DIP,
        #        self.HYPO_DOWN_DIP,
        #        self.DEPTH_TO_TOP))
        #print ("resulting fsx: %f fxy: %f fsz: %s\n" % (self.fsx,
        #                                                self.fsy,
        #                                                self.fsz))

    def calculate_stress(self):
        """
        This function calculates the stress parameters for SDSU based
        on the depth of the fault. These values are calibrated for use
        in Eastern North America
        """
        stress = 16.0 * self.DEPTH_TO_TOP + 225
        stress = stress * 10**6

        return stress

    def __init__(self, a_srcfile=None):
        """
        Set up some parameters for BBToolbox
        """
        self.MAG = None
        self.grid_x = None
        self.grid_y = None
        self.grid_z = 125.0
        self.copy_lf_seismograms = True

        # Parse src file, if given
        if a_srcfile:
            self.parse_src(a_srcfile)

        self.MODALITY = 1
        # GS_FLAG: Don't change it here, override it in the velocity
        # model config file using a 'CODEBASE_SDSU_GS_FLAG = XXX' line
        # 1: Western US (active region),
        # 2: Eastern NA (stable region),
        # 3: Japan
        self.GS_FLAG = 1
        # NGAW_FLAG: Don't change it here, override it in the velocity
        # model config file using a 'CODEBASE_SDSU_NGAW_FLAG = XXX' line
        # 1: NGA-WEST1
        # 2: NGA-WEST2
        self.NGAW_FLAG = 2
        self.KAPPA = 0.04
        self.Q_CODA = 150.0
        self.FDEC = 0.8
        self.AFAC = 41.0
        self.BFAC = 34.0
        self.SOURCE_MECH = "rs"
        self.SOURCE_FUNC = "dreg"
        self.VERBOSE = "on"
        self.TR_SCA = 0.075
        self.STR_FAC = 50.e6

        # 06/10/11: Sandarsh MK
        # Note: Setting FMAX = 20.00 Hz will
        # cause BBtoolbox to produce NaNs in 000 and 090 seismograms.
        self.FMAX = 100.00

        # 09/22/2020: Correlation flag
        self.corr_flag = 2

if __name__ == "__main__":
    BBCODA2 = BBToolboxCfg()
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
