#!/usr/bin/env python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Song RMG Rupture Generator Multi Segment Configuration File
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils

class SongRMGMSCfg(object):
    """
    This class implements the RMG rupture generator configuration options
    """

    def __init__(self, a_srcfiles):
        """
        Sets basic class parameters, then parses a_srcname for more
        information.
        """
        # Available options 'tri', 'rec', 'pliu', 'etinti'
        self.svf_type = 'etinti'
        self.svf_dt = 0.1
        self.num_srcfiles = len(a_srcfiles)
        self.seg_hypocenter = 0
        self.CFGDICT = []
        for a_srcname in a_srcfiles:
            self.CFGDICT.append(bband_utils.parse_src_file(a_srcname))
        # Check SRC files to determine which segment has the hypocenter
        print()
        print("*** Default hypocenter segment: %d" % (self.seg_hypocenter))
        for segment in range(0, self.num_srcfiles):
            if "true_hypo" in self.CFGDICT[segment]:
                if self.CFGDICT[segment]["true_hypo"] == 1:
                    self.seg_hypocenter = segment
                    print("*** Updated hypocenter segment: %d" % (self.seg_hypocenter))
                    break

if __name__ == "__main__":
    RMG_CFG = SongRMGMSCfg(sys.argv[1:])
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
