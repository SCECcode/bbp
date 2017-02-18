#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Getargs class to read arguments and parameters
$Id: syn1D_args.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import sys
from optparse import OptionParser

class syn1DArgs(object):
    """
    argument parser
    """
    def __init__(self):
        usage_string = ("./%ucsyn.py [-s (--simID) -r (--srffile) "
                        "-p (--stationlist) -t (--test) -l (--logfile)]")
        self.parser = OptionParser(usage=usage_string)
        self.getargs()

    def getargs(self):
        self.parser.add_option("-v", "--verbose", action="store_true",
                               help="Set verbose on", dest="verbose_var",
                               default=False)
        self.parser.add_option("-t", "--test", action="store_true",
                               dest="test_mode",
                               help="Specify that commands are printed, but not run.",
                               default=False)
        self.parser.add_option("-l", "--logfile", action="store",
                               type="string", dest="log_file",
                               help="Specify the Log file name",
                               metavar="LOG_FILE_NAME")
        self.parser.add_option("-r", "--srffile", action="store",
                               type="string", dest="a_srf_file",
                               help="An absolute path to the SRF_FILENAME",
                               metavar="SRF_FILENAME")
        self.parser.add_option("-s", "--simID", action="store", type="int",
                               dest="simID",
                               help="Specify the SIMULATION_ID",
                               metavar="SIMULATION_ID")

        (self.options, self.args) = self.parser.parse_args()

        if not self.options.simID:
            self.parser.error('-s --simID is a required parameter')
        #if not self.options.a_srf_file:
        #  self.parser.error('-r --srf is a required parameter')
        #if not self.options.a_station_list:
        #  self.parser.error('-p --stationlist is a required parameter')
        if len(self.args) != 0:
            self.parser.error("incorrect number of arguments")

if __name__ == "__main__":
    GETA = syn1DArgs()
    print(GETA.options)
    sys.exit(0)
