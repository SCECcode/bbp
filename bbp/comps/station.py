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

Class for managing station info in BB Platform
"""
from __future__ import division, print_function

# Import Python modules
import sys

class Station(object):
    """
    Class for storing station information
    """

    lon = None
    lat = None
    scode = None
    vs30 = None
    low_freq_corner = 1.0e-15
    high_freq_corner = 1.0e+15

    # SDSU fields
    x = None
    y = None
    vp = None
    vs = None
    rho = None
    kappa = None

if __name__ == "__main__":
    print("Testing Module: %s" % (sys.argv[0]))
    ME = Station()
    ME.lon = 34.00
    ME.lat = -118.00
    ME.scode = "s003"
    ME.vs30 = 865
    sys.exit(0)
