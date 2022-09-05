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
        self.R_UC_FFSP_EXE = "ffsp_dcf_v2"
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
