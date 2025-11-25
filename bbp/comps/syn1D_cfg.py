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

This module defines the configuration parameters for the UCSB syn_1d program
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
from install_cfg import InstallCfg
import bband_utils
import velocity_models

class Syn1DCfg(object):
    """
    Define the configuration parameters for the syn_1d program
    """

    def __init__(self, vmodel_name, a_srcname=None):

        # Get pointers to all directories
        install = InstallCfg.getInstance()

        # Parse SRC File
        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

        #
        # Name and Path to executables
        self.R_SLL2XY = "statLL2XY"
        self.A_SLL2XY = os.path.join(install.A_UCSB_BIN_DIR,
                                     self.R_SLL2XY)

        self.R_SRF2XY = "srfLL2XYKM"
        self.A_SRF2XY = os.path.join(install.A_UCSB_BIN_DIR,
                                     self.R_SRF2XY)

        self.R_SYN1D = "syn_1d"
        self.A_SYN1D = os.path.join(install.A_UCSB_BIN_DIR,
                                    self.R_SYN1D)

        self.R_CONV = "conv3CompBB"
        self.A_CONV = os.path.join(install.A_UCSB_BIN_DIR,
                                   self.R_CONV)

        self.R_STITCH = "stitch"
        self.A_STITCH = os.path.join(install.A_UCSB_BIN_DIR,
                                     self.R_STITCH)

        #
        # Define name used when input station file is converted into a
        # UC lat/lon version of the station file
        #
        self.R_UC_STATION_FILE = "uc_stations.ll"
        self.R_UC_VS30_FILE = "stations.vs30"
        self.R_UC_SOURCE_MODEL = "source_model.list"
        self.R_FFSP_FILE = "FFSP_OUTPUT.bst"
        self.MAX_STATIONS = 300

        vmodel_obj = velocity_models.get_velocity_model_by_name(vmodel_name)
        if vmodel_obj is None:
            raise IndexError("Cannot find velocity model: %s" %
                             (vmodel_name))

        vmodel_params = vmodel_obj.get_codebase_params('ucsb')
        # Configure needed parameters from velocity model
        if 'LF_VELMODEL' in vmodel_params:
            self.A_UC_LF_VELMODEL = os.path.join(vmodel_obj.base_dir,
                                                 vmodel_params['LF_VELMODEL'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("LF_VELMODEL", vmodel_name))
        if 'HF_VELMODEL' in vmodel_params:
            self.A_UC_HF_VELMODEL = os.path.join(vmodel_obj.base_dir,
                                                 vmodel_params['HF_VELMODEL'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("HF_VELMODEL", vmodel_name))
        if 'GREENBANK' in vmodel_params:
            self.A_UC_GREENBANK = os.path.join(vmodel_obj.base_dir,
                                               vmodel_params['GREENBANK'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("GREENBANK", vmodel_name))
        if 'GREEN_SOIL' in vmodel_params:
            self.A_UC_GREEN_SOIL = os.path.join(vmodel_obj.base_dir,
                                                vmodel_params['GREEN_SOIL'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("GREEN_SOIL", vmodel_name))
        if 'HF_GREENBANK' in vmodel_params:
            self.A_UC_HF_GREENBANK = os.path.join(vmodel_obj.base_dir,
                                                  vmodel_params['HF_GREENBANK'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("HF_GREENBANK", vmodel_name))
        if 'HF_GREEN_SOIL' in vmodel_params:
            self.A_UC_HF_GREEN_SOIL = os.path.join(vmodel_obj.base_dir,
                                                   vmodel_params['HF_GREEN_SOIL'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("HF_GREEN_SOIL", vmodel_name))
        if 'SYN1D_INP_FILE' in vmodel_params:
            self.A_UC_SYN1D_INP_FILE = os.path.join(vmodel_obj.base_dir,
                                                    vmodel_params['SYN1D_INP_FILE'])
        else:
            raise KeyError("%s parameter missing in velocity model %s" %
                           ("SYN1D_INP_FILE", vmodel_name))

if __name__ == "__main__":
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
