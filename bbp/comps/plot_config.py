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

This module contains plot-specific parameters that are shared among
multiple plotting scripts.
"""
from __future__ import division, print_function

dpi = None # *None* will default to the value savefig.dpi in the
           # matplotlibrc file, or use scalar > 0 to specify a
           # different value.

line_width = 1.0 # Line width in points, default is 1.0

plot_seismograms_duration = 100
plot_seismograms_mode = 0  # 0 - Plot "duration" seconds only. Try
                           #     to guess start if seismogram longer
                           #     than "duration" seconds
                           # 1 - Plot x seconds starting at t=0s
                           # 2 - Plot entire seismogram

PLOT_SRF_DEFAULT_CONTOUR_INTERVALS = 2.0 # Contour intervals (in seconds)

# Per-method ranges for PSA and FAS GoF plots
PSA_GOF_FREQ = {'gp': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'sdsu': {'freq_low' : 0.1, 'freq_high' : 50.0},
                'ucsb': {'freq_low' : 0.1, 'freq_high' : 15.0},
                'exsim': {'freq_low' : 0.2, 'freq_high' : 100.0},
                'irikura1': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'irikura2': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'song': {'freq_low' : 0.1, 'freq_high' : 10.0}
                }

FAS_GOF_FREQ = {'gp': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'sdsu': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'ucsb': {'freq_low' : 0.1, 'freq_high' : 9.0},
                'exsim': {'freq_low' : 0.1, 'freq_high' : 40.0},
                'irikura1': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'irikura2': {'freq_low' : 0.1, 'freq_high' : 10.0},
                'song': {'freq_low' : 0.1, 'freq_high' : 10.0}
                }
