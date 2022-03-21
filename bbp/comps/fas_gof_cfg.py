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

This module contains configuration parameters for the FAS GOF generation
"""
from __future__ import division, print_function

# Import Python modules
import sys
import numpy as np

T95 = [ 6.3138, 2.9200, 2.3534, 2.1318, 2.0150, 1.9432, 1.8946,
        1.8595, 1.8331, 1.8125, 1.7959, 1.7823, 1.7709, 1.7613,
        1.7531, 1.7459, 1.7396, 1.7341, 1.7291, 1.7247, 1.7207,
        1.7171, 1.7139, 1.7109, 1.7081, 1.7056, 1.7033, 1.7011,
        1.6991, 1.6973, 1.6955, 1.6939, 1.6924, 1.6909, 1.6896,
        1.6883, 1.6871, 1.6860, 1.6849, 1.6839, 1.6829, 1.6820,
        1.6811, 1.6802, 1.6794, 1.6787, 1.6779, 1.6772, 1.6766,
        1.6759, 1.6753, 1.6747, 1.6741, 1.6736, 1.6730]

def uncer(nstat_read, nval, tmin, tmax, nper, per,
              resid, bias, sigma, sigma0, cl90m, cl90p):
    """
    Calculate GoF parameters
    """
    for i in range(0, nper):
        bias[i] = 0.0
        cl90m[i] = 0.0
        nval[i] = 0

    for i in range(0, nstat_read):
        for j in range(0, nper):
            if per[j] >= tmin[i] and per[j] <= tmax[i]:
                bias[j] = bias[j] + resid[i][j]
                cl90m[j] = cl90m[j] + resid[i][j] * resid[i][j]
                nval[j] = nval[j] + 1

    for j in range(0, nper):
        if nval[j] > 1:
            invn = 1.0 / float(nval[j])
            invn1 = 1.0 / float(nval[j] - 1)

            if nval[j] > 56:
                ttfac = 1.64 * np.sqrt(invn)
            else:
                ttfac = T95[nval[j]-2] * np.sqrt(invn)

            bias[j] = bias[j] * invn
            sigma[j] = np.sqrt(invn1 * (cl90m[j] - nval[j] * bias[j] * bias[j]))
            sigma0[j] = np.sqrt(invn * (cl90m[j]))
            cl90m[j] = bias[j] - sigma[j] * ttfac
            cl90p[j] = bias[j] + sigma[j] * ttfac
        else:
            bias[j] = 0.0
            sigma[j] = 0.0
            sigma0[j] = 0.0
            cl90m[j] = 0.0
            cl90p[j] = 0.0

def resid2uncer_py(residfile, fileroot, comp,
                   nstat, min_cdst, max_cdst):
    """
    Python version of GP resid2uncer code
    """
    min_vs30 = -1e+15
    max_vs30 =  1e+15
    min_xcos = -1e+15
    max_xcos =  1e+15
    min_ycos = -1e+15
    max_ycos =  1e+15
        
    input_file = open(residfile, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        break
    # Get list of periods
    per = line.split()[13:]
    per = [float(period) for period in per]
    nper = len(per)
    # print("periods in file = %d\n" % (nper))

    nstat_read = 0
    tmin = []
    tmax = []
    resid = []
    nval = np.zeros(len(per), dtype=int)
    bias = np.zeros(len(per))
    sigma = np.zeros(len(per))
    sigma0 = np.zeros(len(per))
    cl90m = np.zeros(len(per))
    cl90p = np.zeros(len(per))
        
    # Continue reading input file
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        pieces = line.split()
        vs30 = float(pieces[6])
        cdst = float(pieces[7])
        xcos = float(pieces[8])
        ycos = float(pieces[9])
        rdcomp = pieces[12]

        if (vs30 >= min_vs30 and
            vs30 <= max_vs30 and
            cdst >= min_cdst and
            cdst <= max_cdst and
            xcos >= min_xcos and
            xcos <= max_xcos and
            ycos >= min_ycos and
            ycos <= max_ycos and
            rdcomp == comp):

            # Read this station
            tmin.append(float(pieces[10]))
            tmax.append(float(pieces[11]))

            pieces = pieces[13:]
            pieces = [float(piece) for piece in pieces]
            resid.append(pieces)
                
            nstat_read = nstat_read + 1
            if nstat_read > nstat:
                pass

    input_file.close()

    # Completed reading input file
    # print("nstat_read = %d " % (nstat_read))
    uncer(nstat_read, nval, tmin, tmax, nper, per,
          resid, bias, sigma, sigma0, cl90m, cl90p)

    # Write output to files
    bias_file = open("%s.bias" % (fileroot), 'w')
    for c_per, c_bias in zip(per, bias):
        bias_file.write("%13.5e %13.5e\n" % (c_per, c_bias))
    bias_file.close()

    sigma_file = open("%s.sigma" % (fileroot), 'w')
    for c_per, c_sigma in zip(per, sigma):
        sigma_file.write("%13.5e %13.5e\n" % (c_per, c_sigma))
    sigma_file.close()

    sigma0_file = open("%s.sigma0" % (fileroot), 'w')
    for c_per, c_sigma0 in zip(per, sigma0):
        sigma0_file.write("%13.5e %13.5e\n" % (c_per, c_sigma0))
    sigma0_file.close()

    m90_file = open("%s.m90" % (fileroot), 'w')
    for c_per, c_m90 in zip(per, cl90m):
        m90_file.write("%13.5e %13.5e\n" % (c_per, c_m90))
    m90_file.close()

    p90_file = open("%s.p90" % (fileroot), 'w')
    for c_per, c_p90 in zip(per, cl90p):
        p90_file.write("%13.5e %13.5e\n" % (c_per, c_p90))
    p90_file.close()

class FASGofCfg(object):
    """
    Define the configuration parameters for FAS gof
    """
    def __init__(self):
        self.COMPS_FAS = ["fash1", "fash2", "seas"]
        self.MIN_CDST = 0
        self.MAX_CDST = 25

if __name__ == "__main__":
    MSEIS = FASGofCfg()
    print("Created Test Config Class: %s" % (sys.argv[0]))
