#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Getargs class to read arguments and parameters
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
