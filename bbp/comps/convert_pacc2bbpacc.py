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

Created on Jul 27, 2012
@author: maechlin

Utility to convert PEER format seismogram (acceleration) to SCEC bbp
format.  Expects a pathname to directory containing 3 component PEER
files. Outputs .bbp format file for each set of 3 PEER files. Output
bbp files are based on the station name, which may be a station name
or a record sequence number rsn.
"""
from __future__ import division, print_function

import os
import bbp_formatter

def main():
#
# Logic is
#
# ls the LOMAP dir
# find each RNS
# for each RSN create three names of the three components
# for each RNS create name of output file based on RNS
# call bbp_formatter
#
#
#
    path = './LOMAP'
    listing = os.listdir(path)
    for infile in listing:
        print("current file is: " + infile)

    sta_list = []
    for infile in listing:
        if infile.endswith("_E.acc" or "_N.acc" or "_Z.acc"):
            fields = infile.split("_")
            id = int(fields[0])
            print("Found RECID %d" % (id))
            sta_list.append(id)
        print("next infile")

    for sta in sta_list:
        e_fname = os.path.join(path, "%s_E.acc" % (sta))
        n_fname = os.path.join(path, "%s_N.acc" % (sta))
        z_fname = os.path.join(path, "%s_Z.acc" % (sta))
        bbp_fname = os.path.join(path, "%s.bbp" % (sta))
        print(e_fname)
        print(n_fname)
        print(z_fname)
        print(bbp_fname)
        bbp_formatter.peer2bbp(n_fname, e_fname, z_fname, bbp_fname)
        print("Created BBP file: %s" % (bbp_fname))

if __name__ == '__main__':
    main()
