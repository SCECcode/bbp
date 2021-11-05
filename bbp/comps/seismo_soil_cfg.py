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
