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
