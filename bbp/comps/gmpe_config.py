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

This module contains gmpe-specific configuration parameters that are
used throught the Broadband Platform
"""
from __future__ import division, print_function

# Import Python modules
import numpy

# Import PyNGA modules
import pynga

GMPES = {"NGA-West1" : {"models" : ['AS', 'BA', 'CB', 'CY'],
                        "labels" : ['AS08', 'BA08', 'CB08', 'CY08'],
                        "func" : pynga.NGA08,
                        "periods" : [0.010, 0.011, 0.012, 0.013, 0.015, 0.017,
                                     0.020, 0.022, 0.025, 0.029, 0.032, 0.035,
                                     0.040, 0.045, 0.050, 0.055, 0.060, 0.065,
                                     0.075, 0.085, 0.100, 0.110, 0.120, 0.130,
                                     0.150, 0.170, 0.200, 0.220, 0.240, 0.260,
                                     0.280, 0.300, 0.350, 0.400, 0.450, 0.500,
                                     0.550, 0.600, 0.650, 0.750, 0.850, 1.000,
                                     1.100, 1.200, 1.300, 1.500, 1.700, 2.000,
                                     2.200, 2.400, 2.600, 2.800, 3.000, 3.500,
                                     4.000, 4.400, 5.000, 5.500, 6.000, 6.500,
                                     7.500, 8.500, 10.000]},
         "NGA-West2" : {"models" : ['ASK', 'BSSA', 'CB', 'CY'],
                        "labels" : ['ASK14', 'BSSA14', 'CB14', 'CY14'],
                        "func" : pynga.NGA14,
                        "periods" : [0.010, 0.011, 0.012, 0.013, 0.015, 0.017,
                                     0.020, 0.022, 0.025, 0.029, 0.032, 0.035,
                                     0.040, 0.045, 0.050, 0.055, 0.060, 0.065,
                                     0.075, 0.085, 0.100, 0.110, 0.120, 0.130,
                                     0.150, 0.170, 0.200, 0.220, 0.240, 0.260,
                                     0.280, 0.300, 0.350, 0.400, 0.450, 0.500,
                                     0.550, 0.600, 0.650, 0.750, 0.850, 1.000,
                                     1.100, 1.200, 1.300, 1.500, 1.700, 2.000,
                                     2.200, 2.400, 2.600, 2.800, 3.000, 3.500,
                                     4.000, 4.400, 5.000, 5.500, 6.000, 6.500,
                                     7.500, 8.500, 10.000]},
         "CENA GROUP 1" : {"models" : ['PZT11', 'A0811E', 'S03SCVS'],
                           "labels" : ['PZT11', 'A0811E', 'S03SCVS'],
                           "func" : pynga.CENA1,
                           "periods" : [0.01, 0.04, 0.10, 0.20,
                                        0.40, 1.0, 2.0]}}

def calculate_gmpe(model_group, model_name, Mw, Rjb, Vs30, period, **kw):
    """
    This function is a wrapper for the PyNGA NGA08, NGA14 and CENA1
    functions
    """
    # Figure out which function to call
    gmpe_group = GMPES[model_group]
    func = gmpe_group['func']
    # CENA GROUP 1 needs Rrup
    if model_group == "CENA GROUP 1":
        Rrup = kw['Rrup']
        median = func(model_name, Mw, Rjb, Rrup, period)
    else:
        median, _, _, _ = func(model_name, Mw, Rjb, Vs30, period, **kw)
        median = (median.tolist())[0]

    # Return median
    return median

def average_gmpe(station, a_src_gmpe_file, a_avg_out_file):
    """
    This function reads the input gmpe file and averages all gmpes
    so that the output file can be plotted in maps/distance plots.
    """
    input_file = open(a_src_gmpe_file, 'r')
    output_file = open(a_avg_out_file, 'w')
    output_file.write("#station: %s\n" % (station))
    output_file.write("#period avg_gmpe avg_gmpe avg_gmpe\n")
    for line in input_file:
        if line.startswith("#") or line.startswith("%"):
            # Skip comments
            continue
        line = line.strip()
        tokens = line.split()
        # Convert to float
        tokens = [float(token) for token in tokens]
        # Remove period from list
        period = tokens.pop(0)
        # Calculate average
        gmpe_avg = numpy.mean(tokens)
        # Write output
        output_file.write("%10.4f    %10.5e    %10.5e    %10.5e\n" %
                          (period, gmpe_avg, gmpe_avg, gmpe_avg))
    # Close files
    input_file.close()
    output_file.close()
