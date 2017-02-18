"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Created on Jul 27, 2012
@author: maechlin

Utility to convert PEER format seismogram (acceleration) to SCEC bbp
format.  Expects a pathname to directory containing 3 component PEER
files. Outputs .bbp format file for each set of 3 PEER files. Output
bbp files are based on the station name, which may be a station name
or a record sequence number rsn.

$Id: convert_pacc2bbpacc.py 1730 2016-09-06 20:26:43Z fsilva $
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
