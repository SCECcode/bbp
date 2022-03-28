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

This config class will encapsulate the configuration parameters
needed to run a simulation. Programs can derive specific
configuration sets from this base class to suppor their own programs.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

class MatchCfg(object):
    """
    Define the configuration parameters for the Match program
    """
    def __init__(self):
        self.HF_FHI = 1.0
        self.HF_FLO = "1.0e+15"
        self.HF_ORD = 4
        self.HF_TSTART = 0.0

        self.LF_FHI = 0.0
        self.LF_FLO = 1.0
        self.LF_ORD = 4
        self.LF_TSTART = 0.0

        self.PHASE = 0
        self.MATCH_METHOD = 2

        self.COMPS = ['000', '090', 'ver']

        self.FILTLIST = "filtmatchlist1"

if __name__ == "__main__":
    print("Test Config Class: %s" % os.path.basename((sys.argv[0])))
