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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
# import random

# Import Broadband modules
import bband_utils

class IrikuraGenSrfCfg(object):
    """
    Define the configuration parameters for the Irikura rupture generator
    """

    def __init__(self, a_srcname=None):
        """
        Sets basic class parameters, then parses a_srcname for more information.
        """
        self.VS = 3.5
        self.DT = 0.025
        self.DENS = 2.7
        self.GENSRF = "gen_srf"

        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

if __name__ == "__main__":
    IRIKURA_GEN_SRF_CFG = IrikuraGenSrfCfg()
    print("Created Test Config Class: %s" % os.path.basename((sys.argv[0])))
