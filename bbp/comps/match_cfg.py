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

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

class MatchCfg(object):
    """
    Define the configuration parameters for the Match program
    """
    def __init__(self):
        self.HF_FHI = 1.0
        self.HF_FLO = "1.0e+15"
        self.HF_ORD = 4
        self.HF_TSTART = 0.0

        self.LF_FHI = 0.0
        self.LF_FLO = 1.0
        self.LF_ORD = 4
        self.LF_TSTART = 0.0

        self.PHASE = 0
        self.MATCH_METHOD = 2

        self.COMPS = ['000', '090', 'ver']

        self.FILTLIST = "filtmatchlist1"

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename((sys.argv[0])))
