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

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import velocity_models

class JbsimCfg(object):
    """
    Define the configuration parameters for the Jbrun program
    """
    def __init__(self, vmodel_name):

        vmodel_obj = velocity_models.get_velocity_model_by_name(vmodel_name)
        if vmodel_obj is None:
            raise IndexError("Cannot find velocity model: %s" %
                             (vmodel_name))

        vmodel_params = vmodel_obj.get_codebase_params('gp')
        # Configure needed parameters from velocity model
        if 'GF_DIR' in vmodel_params:
            self.A_GP_GF_DIR = os.path.join(vmodel_obj.base_dir,
                                            vmodel_params['GF_DIR'])
        else:
            raise KeyError("GF_DIR parameter missing in velocity model %s" %
                           (vmodel_name))
        if 'GF_NAME' in vmodel_params:
            self.GF_NAME = vmodel_params['GF_NAME']
        else:
            raise KeyError("GF_NAME parameter missing in velocity model %s" %
                           (vmodel_name))
        if 'GF_LOCS' in vmodel_params:
            self.GF_LOCS = os.path.join(vmodel_obj.base_dir,
                                        vmodel_params['GF_LOCS'])
        else:
            raise KeyError("GF_LOCS parameter missing in velocity model %s" %
                           (vmodel_name))
        if 'GF_TIMES' in vmodel_params:
            self.GF_TIMES = os.path.join(vmodel_obj.base_dir,
                                         vmodel_params['GF_TIMES'])
        else:
            raise KeyError("GF_TIMES parameter missing in velocity model %s" %
                           (vmodel_name))
        if 'GF_DT' in vmodel_params:
            self.MIN_GFDT = float(vmodel_params['GF_DT'])
        else:
            raise KeyError("GF_DT parameter missing in velocity model %s" %
                           (vmodel_name))
        # Now, look for optional parameters in the velocity model
        if 'MAX_GFNT' in vmodel_params:
            self.MAX_GFNT = int(vmodel_params['MAX_GFNT'])
        else:
            self.MAX_GFNT = 4096
        if 'NTOUT' in vmodel_params:
            self.NTOUT = int(vmodel_params['NTOUT'])
        else:
            self.NTOUT = 4096

        # Now, configure other paramteres
        self.GF_SWAP_BYTES = 0
        self.DTOUT = self.MIN_GFDT
        self.COMPS = ["000", "090", "ver"]

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename((sys.argv[0])))
