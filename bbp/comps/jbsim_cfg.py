#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
$Id: jbsim_cfg.py 1722 2016-08-26 21:50:46Z fsilva $
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
        self.VREF = 865

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename((sys.argv[0])))
