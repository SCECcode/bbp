#!/usr/bin/env python
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

Class to create sequence number generator. Uses a small file text, but will
work with defaults if it has to.
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
