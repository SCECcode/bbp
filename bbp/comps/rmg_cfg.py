#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

RMG Rupture Generator Configuration File
$Id: rmg_cfg.py 1762 2016-09-20 21:14:57Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils

class RMGCfg(object):
    """
    This class implements the RMG rupture generator configuration options
    """

    def __init__(self, a_srcname=None):
        """
        Sets basic class parameters, then parses a_srcname for more
        information.
        """
        # Available options 'tri', 'rec', 'pliu', 'etinti'
        self.svf_type = 'etinti'
        self.svf_dt = 0.1

        if a_srcname:
            self.CFGDICT = bband_utils.parse_src_file(a_srcname)

if __name__ == "__main__":
    RMG_CFG = RMGCfg()
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
