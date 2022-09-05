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

Configuration file for SeismoSoil site response module
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import velocity_models

class SeismoSoilCfg(object):
    """
    Define the configuration parameters for SeismoSoil
    """
    def __init__(self, vmodel_name, method):

        vmodel_obj = velocity_models.get_velocity_model_by_name(vmodel_name)
        if vmodel_obj is None:
            raise IndexError("Cannot find velocity model: %s" %
                             (vmodel_name))

        if method.lower() == "ucsb":
            vmodel_params = vmodel_obj.get_codebase_params('ucsb')
        elif method.lower() == "exsim":
            vmodel_params = vmodel_obj.get_codebase_params('exsim')
        elif method.lower() == "sdsu":
            vmodel_params = vmodel_obj.get_codebase_params('sdsu')
        else:
            # For now...
            vmodel_params = vmodel_obj.get_codebase_params('gp')

        # Read reference velocities for LF and HF components, use defaults
        # values if not found so that the code will still work without GP GFs
        if 'LF_VREF' in vmodel_params:
            self.LF_VREF = int(vmodel_params['LF_VREF'])
        else:
            raise ValueError("Cannot find LF_VREF for method %s!" % (method))
        if 'HF_VREF' in vmodel_params:
            self.HF_VREF = int(vmodel_params['HF_VREF'])
        else:
            raise ValueError("Cannot find HF_VREF for method %s!" % (method))

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename(sys.argv[0]))
