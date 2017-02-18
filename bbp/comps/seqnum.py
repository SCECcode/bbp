#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Class to create sequence number generator. Uses a small file text, but will
work with defaults if it has to.
$Id: seqnum.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import sys
import time
import datetime

def get_seq_num():
    """
    Simple class for creating sequence numbers
    Truncate epoch time to 7 digits which is about one month
    """
    t = datetime.datetime.now()
    mt = time.mktime(t.timetuple())
    nextnum = int(mt)
    retval = nextnum % 10000000
    return retval

if __name__ == "__main__":
    val = get_seq_num()
    val2 = get_seq_num()
    print("found : %d and %d" % (val, val2))
    sys.exit(0)
