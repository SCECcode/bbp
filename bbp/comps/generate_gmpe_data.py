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

Workflow using pynga for broadband platform
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import pynga and its utilities
from pynga import *
from pynga.utils import *

if len(sys.argv) < 3:
    print("Usage %s: src_file stat_file" % (sys.argv[0]))
    sys.exit(0)

source_file = sys.argv[1]
stafile = sys.argv[2]

# ==============================
# 1. Read Event Info from SRC file
# ==============================
cfg_dict = {}
src_file = open(source_file, 'r')
for line in src_file:
    # Read file, line by line
    line = line.strip()
    if line.startswith('#'):
        # Skip comments
        continue
    ml = line.split('=')
    key = ml[0].strip().upper()
    val = ml[1].strip()
    cfg_dict[key] = float(val)
src_file.close()

# Look for the keys that we need
M = cfg_dict['MAGNITUDE']
L = cfg_dict['FAULT_LENGTH']
dl = cfg_dict['DLEN']
W = cfg_dict['FAULT_WIDTH']
dw = cfg_dict['DWID']
ztor = cfg_dict['DEPTH_TO_TOP']
strike = cfg_dict['STRIKE']
rake = cfg_dict['RAKE']
dip = cfg_dict['DIP']
lat0 = cfg_dict['LAT_TOP_CENTER']
lon0 = cfg_dict['LON_TOP_CENTER']
hypoAS = cfg_dict['HYPO_ALONG_STK']
hypoDD = cfg_dict['HYPO_DOWN_DIP']

# information needed to get distance parameters
origin = lon0, lat0    # center
Dims = L, dl, W, dw, ztor
Mech = strike, dip, rake

# ======================
# 2. Read periods list
# ======================

# Periods in seconds (could be any values in between 0.01 to 10)
periods = [0.010, 0.011, 0.012, 0.013, 0.015, 0.017, 0.020, 0.022,
           0.025, 0.029, 0.032, 0.035, 0.040, 0.045, 0.050, 0.055,
           0.060, 0.065, 0.075, 0.085, 0.100, 0.110, 0.120, 0.130,
           0.150, 0.170, 0.200, 0.220, 0.240, 0.260, 0.280, 0.300,
           0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.750,
           0.850, 1.000, 1.100, 1.200, 1.300, 1.500, 1.700, 2.000,
           2.200, 2.400, 2.600, 2.800, 3.000, 3.500, 4.000, 4.400,
           5.000, 5.500, 6.000, 6.500, 7.500, 8.500, 10.000]

# For nga_model, you can use 'CB','CY', and 'AS' as well, or loop over
# all four NGA models
nga_models = ['AS', 'BA', 'CB', 'CY']

# =====================
# 3. Read Station list and compute distance parameters
# Vs30, Z1.0 and Z2.5 are also needed from Station List
# =====================
stations = open(stafile)
for line in stations:
    line = line.strip()
    # Skip commments
    if line.startswith('#'):
        continue
    tokens = line.split()
    if len(tokens) < 5:
        # Not enough items
        continue

    station = tokens[2]

    # Station location
    SiteGeom = [float(tokens[0]), float(tokens[1]), 0.0]
    FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, dummy1, dummy2 = FaultTraceGen( origin, Dims, Mech )
    Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom, FaultTrace1,
                                                 UpperSeisDepth,
                                                 LowerSeisDepth,
                                                 AveDip) # in km

    # assuming you get the following site condition from station list:
    # just for one, you could use a list of station parameters
    Vs30 = float(tokens[4]) # in m/s
    #Z10 = 1000    # in meter
    #Z25 = 3.0     # in km
    Z10 = None    # Let Pynga calculate it
    Z25 = None    # Let Pynga calculate it

    # =====================
    # 4. Compute PSA for the given station
    # =====================

    # For a given event and a given station, this list has
    # period-dependent PSA values
    StationMedian = []
    for ip in xrange(len(periods)):
        period = periods[ip]
        period_medians = []
        for nga_model in nga_models:
            # this input list should be enough based on what BBP has
            median, sigmaT, tau, sigma = NGA08(nga_model, M, Rjb, Vs30,
                                               period,
                                               rake=rake, dip=dip, W=W,
                                               Ztor=ztor, Rrup=Rrup, Rx=Rx,
                                               Z10=Z10, Z25=Z25)
            # median in (g) for PSA and PGA
            period_medians.append(median)
        StationMedian.append((period, period_medians))

    # =====================
    # 5. Output results to file
    # =====================

    outfile = open("%s-gmpe.ri50" % (station), 'w')
    outfile.write("#station: %s\n" % (station))
    outfile.write("#period AS BA CB CY\n")
    for item in StationMedian:
        period = item[0]
        vals = item[1]
        out_str = "%.4f" % (period)
        for method in vals:
            out_str = out_str + "\t%.6f" % (method)
        outfile.write("%s\n" % (out_str))
    outfile.close()

# ======================
# 6. Write metadata file
# ======================
args = " ".join([os.path.abspath(arg) for arg in sys.argv])
metafile = open("params.txt", 'w')
metafile.write("#GMPEworkflow.py inputs\n")
metafile.write("$ %s\n\n" % (args))
metafile.close()
