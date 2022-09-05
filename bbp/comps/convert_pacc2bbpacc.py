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
