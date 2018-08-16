#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module contains configuration parameters for the GOF generation
"""
from __future__ import division, print_function

# Import Python modules
import sys

class GPGofCfg(object):
    """
    Define the configuration parameters for gof
    """
    def __init__(self):
        self.COMPS = ["090", "000", "avgh"]
        self.COMPS_PSA5 = ["psa5n", "psa5e", "rotd50"]
        self.MIN_CDST = 0
        self.MAX_CDST = 25

if __name__ == "__main__":
    MSEIS = GPGofCfg()
    print("Created Test Config Class: %s" % (sys.argv[0]))
