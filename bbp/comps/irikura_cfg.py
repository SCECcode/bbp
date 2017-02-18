#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
# import random

# Import Broadband modules
import bband_utils

class IrikuraCfg(object):
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
    IRIKURA_CFG = IrikuraCfg()
    print("Created Test Config Class: %s" % os.path.basename((sys.argv[0])))
