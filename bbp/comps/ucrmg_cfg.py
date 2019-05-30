#!/usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module is used to configure the UCSB rupture generator parameters.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
from install_cfg import InstallCfg
import bband_utils
import velocity_models

class UCrmgCfg(object):
    """
    Define the configuration parameters for the UC SB Rupture Model Generator
    """

    def __init__(self, vmodel_name, a_srcname=None):
        install = InstallCfg.getInstance()
        #
        # Name and Path to executable
        #
        self.R_UC_FFSP_EXE = "ffsp_v2"
        self.A_UC_FFSP_EXE = os.path.join(install.A_UCSB_BIN_DIR,
                                          self.R_UC_FFSP_EXE)
        self.FFSP_OUTPUT_PREFIX = "FFSP_OUTPUT"

        self.FMAX = 50.0 # Nyquist -- use 50 for 100Hz

        vmodel_obj = velocity_models.get_velocity_model_by_name(vmodel_name)
        if vmodel_obj is None:
            raise IndexError("Cannot find velocity model: %s" %
                             (vmodel_name))

        vmodel_params = vmodel_obj.get_codebase_params('ucsb')
        # Configure DT based on information from velocity model
        if 'GF_DT' in vmodel_params:
            self.DT = float(vmodel_params['GF_DT'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("GF_DT", vmodel_name))

        # Other region-specific parameters
        if 'RV_AVG' in vmodel_params:
            self.RV_AVG = float(vmodel_params['RV_AVG'])
        else:
            self.RV_AVG = 2.5

        if 'TP_TR' in vmodel_params:
            self.TP_TR = float(vmodel_params['TP_TR'])
        else:
            self.TP_TR = 0.1

        if 'LF_VELMODEL' in vmodel_params:
            self.A_UC_LF_VELMODEL = os.path.join(vmodel_obj.base_dir,
                                                 vmodel_params['LF_VELMODEL'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("LF_VELMODEL", vmodel_name))

        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

            # RV_AVG is optional!
            # If SRC file has it, it overrides the region and the default values
            if "rv_avg" in self.CFGDICT:
                self.RV_AVG = self.CFGDICT["rv_avg"]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please specify velocity model name!")
        sys.exit(1)
    CFG = UCrmgCfg(sys.argv[1])
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
