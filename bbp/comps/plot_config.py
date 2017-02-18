#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This module contains plot-specific parameters that are shared among
multiple plotting scripts.
$Id: plot_config.py 1785 2017-01-24 19:40:10Z fsilva $
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
