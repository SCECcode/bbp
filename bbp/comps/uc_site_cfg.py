#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Configuration file for the UCSB site response module
"""

from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg

class UCSiteCfg(object):
    """
    Define the configuration parameters for the UCSB site response program
    """

    def __init__(self, a_srcname=None):

        # Get pointers to all directories
        install = InstallCfg.getInstance()

        # Parse SRC File
        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

        #
        # Name and Path to executables
        self.R_UC_DECON_EXE = "deconvBBP"
        self.A_UC_DECON_EXE = os.path.join(install.A_UCSB_BIN_DIR,
                                           self.R_UC_DECON_EXE)

        self.R_SLL2XY = "statLL2XY"
        self.A_SLL2XY = os.path.join(install.A_UCSB_BIN_DIR,
                                     self.R_SLL2XY)

        self.R_STITCH = "stitchBBP"
        self.A_STITCH = os.path.join(install.A_UCSB_BIN_DIR,
                                     self.R_STITCH)

        #
        # Define name used when input station file is converted into a UC lat/lon version
        # of the station file
        #
        self.R_UC_STATION_FILE = "uc_stations.ll"
        self.R_UC_VS30_FILE = "stations.vs30"
        self.COMPS = ['000', '090', 'ver']

if __name__ == "__main__":
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
