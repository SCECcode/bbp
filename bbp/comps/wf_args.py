#!/bin/env python
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
