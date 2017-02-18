#!/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Getargs class to read arguments and parameters
"""
from __future__ import division, print_function

# Import Python modules
import sys
from optparse import OptionParser

class wf_args(object):
    """
    Workflow arguement parser. Either no parameters are specified, in which case all defaults are used,
    or else all parameters must be specified"
    """

    def __init__(self):
        usage_string = "./%wf_xxx_xxx_xxx.py [-i --simID -s (--srcfile) -v (--velocitymodel) -l (--stationlist) ]"
        self.parser = OptionParser(usage=usage_string)
        self.parser.add_option("-i", "--simID", action="store", type="int",
                               dest="simID",
                               help="Specify the Simulation ID", metavar="SIM_ID")
        self.parser.add_option("-s", "--srcfile", action="store", type="string",
                               dest="src_file",
                               help="Specify the Source file name",
                               metavar="SOURCE_FILE_NAME")
        self.parser.add_option("-v", "--velmodel", action="store",
                               type="string", dest="velmodel",
                               help="Specify the Velocity Model file name",
                               metavar="VELMODEL_FILE_NAME")
        self.parser.add_option("-l", "--stationlist", action="store",
                               type="string", dest="station_file",
                               help="Specify the Station List file name",
                               metavar="STATION_FILE_NAME")
        (self.options, self.args) = self.parser.parse_args()

        if not self.options.simID:
            self.parser.error('-i --simID is a required parameter')
        if not self.options.src_file:
            self.parser.error('-s --srcfile is a required parameter')
        if not self.options.velmodel:
            self.parser.error('-v --velmodel is a required parameter')
        if not self.options.station_file:
            self.parser.error('-l --stationfile is a required parameter')

if __name__ == "__main__":
    GETA = wf_args()
    print(GETA.options)
    sys.exit(0)
