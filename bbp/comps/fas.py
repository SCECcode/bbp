#!/usr/bin/env python3
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

This module calculates FAS for both simulated and observed timeseries.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import time
import shutil
import matplotlib as mpl
mpl.use('AGG')
import pylab
import numpy as np

# Import Broadband modules
import bband_utils
import install_cfg
import plot_config
from station_list import StationList

# Import plot config file
import plot_config

COMPS = ['000', '090', 'ver']

def create_boore_asc2smc(control_file, input_file,
                         data_column, num_headers,
                         extension_string):
    """
    This function creates the control file for the asc2smc converter tool
    """
    ctl_file = open(control_file, 'w')
    ctl_file.write("!Control file for ASC2SMC      ! first line\n")
    ctl_file.write("! Revision of program involving a change in the "
                   "control file on this date:\n")
    ctl_file.write("   02/02/12\n")
    ctl_file.write("!Name of summary file:\n")
    ctl_file.write(" asc2smc.sum\n")
    ctl_file.write("!n2skip (-1=headers preceded by !; 0=no headers; "
                   "otherwise number of headers to skip)\n")
    ctl_file.write(" %d\n" % (num_headers))
    ctl_file.write("!write headers to smc file "
                   "(even if n2skip > 0)? (Y/N)\n")
    ctl_file.write(" Y\n")
    ctl_file.write("!sps (0.0 = obtain from input file)\n")
    ctl_file.write(" 0\n")
    ctl_file.write("!N columns to read, column number for "
                   "time and data columns \n")
    ctl_file.write("!  (for files made using blpadflt, period is in "
                   "column 1 and sd, pv, pa, rv, \n")
    ctl_file.write("!  aa are in columns 2, 3, 4, 5, 6, respectively)\n")
    ctl_file.write("! Note: if sps .ne. 0.0, then column number for time "
                   "is ignored (but a placeholder is\n")
    ctl_file.write("! still needed--e.g., 1 1 1 (read one column, which "
                   "contains the data; 1 20 1 would be the same)\n")
    ctl_file.write("! But note: if the data are not in the first column, "
                   "but only the data column is to be read\n")
    ctl_file.write("! (because sps will be used to establish "
                   "the time values),\n")
    ctl_file.write("! then ncolumns must be the column corresponding to "
                   "the data.  For example, assume that\n")
    ctl_file.write("! the data are in column 3 and that columns 1 and 2 "
                   "contain time and some other variable, but\n")
    ctl_file.write("! the time column is not to be used (perhaps because "
                   "accumulated error in creating the column\n")
    ctl_file.write("! leads to a slight shift in the time values).  "
                   "Then the input line should be:\n")
    ctl_file.write("!  3 1 3\n")
    ctl_file.write("!\n")
    ctl_file.write("! This program assumes one data point per row; if "
                   "there are more points (as, for example,\n")
    ctl_file.write("! in files with N points per line), "
                   "use the program wrapped2asc).\n")
    ctl_file.write("!\n")
    ctl_file.write(" 3 1 %d\n" % (data_column))
    ctl_file.write("!Xfactr\n")
    ctl_file.write(" 1.0\n")
    ctl_file.write("!Read input format (used if the format is such that "
                   "the values are not separated by spaces,\n")
    ctl_file.write("!in which case a free format cannot be "
                   "used for input)?\n")
    ctl_file.write("  N\n")
    ctl_file.write("!If yes, specify a format; if not, "
                   "still need a placeholder\n")
    ctl_file.write(" (3e13.5)\n")
    ctl_file.write("!For output, use old (standard) smc format or new\n")
    ctl_file.write('!higher precision format.   Specify "high" for\n')
    ctl_file.write("!high precision; any other word defaults to standard\n")
    ctl_file.write("!precision (but some word is needed as "
                   "a placeholder, even if\n")
    ctl_file.write("!standard precision is desired).\n")
    ctl_file.write(" high\n")
    ctl_file.write("!String to append to input file name "
                   "for the output filename.\n")
    ctl_file.write(" %s\n" % (extension_string))
    ctl_file.write('!Input file name (time,data pairs; "stop" in any '
                   'column to quit):\n')
    ctl_file.write("%s\n" % (input_file))
    ctl_file.write("STOP\n")
    ctl_file.close()

def create_boore_smc2fs2(control_file, input_file, name_string):
    """
    This function creates the control file for the smc2fs2 FAS tool
    """
    ctl_file = open(control_file, 'w')
    ctl_file.write('!Control file for program SMC2FS2\n')
    ctl_file.write('! Revision of program involving a change in the control '
                   'file on this date:\n')
    ctl_file.write('   03/10/10\n')
    ctl_file.write('! As many comment lines as desired, each '
                   'starting with "!"\n')
    ctl_file.write('! The string "pp:" indicates a new set '
                   'of processing parameters\n')
    ctl_file.write('! to be applied to the following smc files.  '
                   'The parameters are given on the\n')
    ctl_file.write('! lines following "pp:", until the next "pp:" line '
                   'or until "stop" is \n')
    ctl_file.write('! encountered.\n')
    ctl_file.write('! NOTE: Use the tapers with caution, '
                   'choosing them so that important signal\n')
    ctl_file.write('! is not reduced by the tapering.  '
                   'This can be particularly a problem with \n')
    ctl_file.write('! analog data from relatively small earthquakes '
                   'that triggered near the \n')
    ctl_file.write('! S-wave arrival.  \n')
    ctl_file.write('!\n')
    ctl_file.write('! -----------------------------------------'
                   '------------------------------------\n')
    ctl_file.write('!\n')
    ctl_file.write('! Meaning of smoothing input parameters\n')
    ctl_file.write('!\n')
    ctl_file.write('! NO SMOOTHING\n')
    ctl_file.write('! itype = 0\n')
    ctl_file.write('! SMOOTHING OVER EQUALLY SPACED FREQUENCIES\n')
    ctl_file.write('! itype = 1: box weighting function\n')
    ctl_file.write('!   smooth_param = width of box weighting function (Hz)\n')
    ctl_file.write('! itype = 2: triangular weighting function\n')
    ctl_file.write('!   smooth_param = width of triangular '
                   'weighting function (Hz)\n')
    ctl_file.write('! SMOOTHING OVER LOGARITHMICALLY SPACED FREQUENCIES\n')
    ctl_file.write('! itype = 3: box weighting function\n')
    ctl_file.write('!   smooth_param = xi, which is the fraction of '
                   'a decade for the\n')
    ctl_file.write('!                  box weighting function \n')
    ctl_file.write('! itype = 4: triangular weighting function\n')
    ctl_file.write('!   smooth_param = xi, which is the fraction of '
                   'a decade for the\n')
    ctl_file.write('!                  triangular weighting function \n')
    ctl_file.write('! itype = 5: Konno and Ohmachi weighting function '
                   '(see BSSA 88, 228-241)\n')
    ctl_file.write('!   smooth_param = xi, which is the fraction '
                   'of a decade for which\n')
    ctl_file.write('!                  the Konno and Ohmachi weighting '
                   'function is greater\n')
    ctl_file.write('!                  than 0.043.(it is related to\n')
    ctl_file.write('!                  their smoothing parameter b '
                   'by the equation\n')
    ctl_file.write('!                  b = 4.0/smooth_param, so we have '
                   'this correspondence between\n')
    ctl_file.write('!                  b and smooth_param\n')
    ctl_file.write('!                         b smooth_param \n')
    ctl_file.write('!                        10         0.40\n')
    ctl_file.write('!                        20         0.20\n')
    ctl_file.write('!                        40         0.10\n')
    ctl_file.write('!                  \n')
    ctl_file.write('!                  b = 40 seems to be commonly used, '
                   'but I do not think that it\n')
    ctl_file.write('!                  gives enough smoothing; '
                   'I PREFER SMOOTH_PARAM = 0.2, \n')
    ctl_file.write('!                  corresponding to b = 20. \n')
    ctl_file.write('!\n')
    ctl_file.write('! ipow = power of FAS to be smoothed '
                   '(2 = smoothing energy spectrum)\n')
    ctl_file.write('!\n')
    ctl_file.write('! df_smooth: Note: need df_smooth for '
                   'linearly-spaced smoothers, \n')
    ctl_file.write('! and generally it should be the df from the fft.  '
                   'For general x data, it is\n')
    ctl_file.write('! the spacing between x values, assumed to be constant,  '
                   'The reason for\n')
    ctl_file.write('! including it as an input parameter is to "fool" the\n')
    ctl_file.write('! program to do smoothing over a specified '
                   'number of points by\n')
    ctl_file.write('! setting df_smooth = 1 and smooth_param = number '
                   'of points (including \n')
    ctl_file.write('! points with zero weight at ends; e.g., '
                   'smooth_param = 5 will \n')
    ctl_file.write('! give a smoother with weights 0, 1/4, 2/4, 1/4, 0; '
                   'smooth_param\n')
    ctl_file.write('! should be odd).\n')
    ctl_file.write('!\n')
    ctl_file.write('! ------------------------------------'
                   '-----------------------------------------\n')
    ctl_file.write('! Meaning of frequency specification parameters:\n')
    ctl_file.write('!\n')
    ctl_file.write('!SPECIFY_FREQUENCIES? (y/n):\n')
    ctl_file.write('! <enter Y or N>\n')
    ctl_file.write('!FREQUENCY SPECIFICATION: \n')
    ctl_file.write('!  If specify_frequencies = Y, then enter the \n')
    ctl_file.write('!    number of frequencies, freq(1), freq(2)..., '
                   'freq(nfreq)\n')
    ctl_file.write('!  If specify_frequencies = N, then enter \n')
    ctl_file.write('!    f_low, f_high, log-spaced (0=N, 1=Y), freq_param\n')
    ctl_file.write('!         if freq_param = 0.0, there is no interpolation, '
                   'and the FFT frequencies \n')
    ctl_file.write('!            are used between f_low and f_high '
                   '(log-spaced is ignored).\n')
    ctl_file.write('!         if freq_param /= 0.0 and log-spaced = 0, '
                   'then freq_param is the spacing of the\n')
    ctl_file.write('!            interpolated frequencies '
                   'between f_low and f_high\n')
    ctl_file.write('!         if freq_param /= 0.0 and log-spaced = 1, '
                   'then freq_param is the number of \n')
    ctl_file.write('!            interpolated frequencies between f_low and '
                   'f_high (NOTE: f_low must be > 0.0)\n')
    ctl_file.write('! ---------------------------------------'
                   '--------------------------------------\n')
    ctl_file.write('!\n')
    ctl_file.write('!Name of summary file:\n')
    ctl_file.write(' smc2fs2.sum\n')
    ctl_file.write('PP: new set of parameters\n')
    ctl_file.write('!tskip, tlength\n')
    ctl_file.write('   0.0 2000.0\n')
    ctl_file.write('!dc_remove?\n')
    ctl_file.write('  .true.        \n')
    ctl_file.write('!Length of taper at beginning and end of time series, '
                   'before adding zeros\n')
    ctl_file.write('! to make the number of points in '
                   'the record a power of two.\n')
    ctl_file.write(' 0.0 0.0\n')
    ctl_file.write('!signnpw2(<0, backup for npw2, no zpad):\n')
    ctl_file.write(' +1.0\n')
    ctl_file.write('!smoothing: itype, ipow, df_smooth '
                   '(0 = FFT df), smooth_param\n')
    ctl_file.write('! (see above for the meaning of these input parameters):\n')
    ctl_file.write('   0 1 0.0 0.20\n')
    ctl_file.write('!SPECIFY_FREQUENCIES? (y/n):\n')
    ctl_file.write('  N\n')
    ctl_file.write('!FREQUENCY SPECIFICATION\n')
    ctl_file.write('   0.01 100.0 0 0.0 \n')
    ctl_file.write('!character string to append to filename:\n')
    ctl_file.write('   %s\n' % (name_string))
    ctl_file.write('!Output in smc format (Y,N)?\n')
    ctl_file.write('! ***IMPORTANT NOTE: Output cannot be in smc '
                   'format if use log-spaced \n')
    ctl_file.write('! frequencies because programs such as smc2asc '
                   'have not been modified\n')
    ctl_file.write('! to deal with log-spaced frequency.\n')
    ctl_file.write(' n\n')
    ctl_file.write('!Files to process:\n')
    ctl_file.write('%s\n' % (input_file))
    ctl_file.write('stop\n')
    ctl_file.close()

def add_extra_points(input_bbp_file, output_bbp_file, num_points):
    """
    Add num_points data points at the end of the input_bbp_File,
    writing the result to output_bbp_file
    """
    bbp_dt = None
    bbp_t1 = None
    bbp_t2 = None
    input_file = open(input_bbp_file, 'r')
    output_file = open(output_bbp_file, 'w')

    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            output_file.write("%s\n" % (line))
            continue
        pieces = line.split()
        cur_dt = float(pieces[0])
        if bbp_dt is None:
            if bbp_t1 is None:
                bbp_t1 = cur_dt
            elif bbp_t2 is None:
                bbp_t2 = cur_dt
                bbp_dt = bbp_t2 - bbp_t1
        # Still write the line to output file
        output_file.write("%s\n" % (line))

    # Close input file
    input_file.close()

    if bbp_dt is None:
        raise bband_utils.ParameterError("Cannot find DT in %s!" %
                                         (input_bbp_file))
    for _ in range(0, num_points):
        cur_dt = cur_dt + bbp_dt
        output_file.write("%5.7f   %5.9e   %5.9e    %5.9e\n" %
                          (cur_dt, 0.0, 0.0, 0.0))

    output_file.close()

def read_bbp_samples(bbp_file):
    """
    Reads BBP file and return number of data points in the timeseries
    """
    number_of_samples = 0

    input_file = open(bbp_file)
    for line in input_file:
        line = line.strip()
        if line.startswith("#") or line.startswith("%"):
            continue
        # Count this data point
        number_of_samples = number_of_samples + 1
    input_file.close()

    return number_of_samples

def read_bbp_dt(bbp_file):
    """
    Reads BBP file and returns the DT
    """
    # Pick DT from these files
    bbp_dt = None
    input_file = open(bbp_file)
    for line in input_file:
        line = line.strip()
        if line.startswith("#") or line.startswith("%"):
            continue
        # Got to first timestamp. Now, pick two consecutive
        # timestamps values
        bbp_t1 = float(line.strip().split()[0])
        bbp_t2 = float(next(input_file).strip().split()[0])
        # Subtract the two times
        bbp_dt = bbp_t2 - bbp_t1
        # All done!
        break
    input_file.close()

    if bbp_dt is None:
        raise bband_utils.ParameterError("Cannot find DT in %s!" %
                                         (bbp_file))

    return bbp_dt

def read_fas_file(fas_file):
    """
    Reads FAS file and returns freq and fas arrays
    """
    freqs = []
    fas = []

    # Read input file
    input_file = open(fas_file, 'r')
    # Skip headers
    for line in input_file:
        line = line.strip()
        # skip blank lines
        if not line:
            continue
        if line.startswith("freq"):
            break
    for line in input_file:
        line = line.strip()
        # skip blank lines
        if not line:
            continue
        pieces = line.split()
        pieces = [float(piece) for piece in pieces]
        freqs.append(pieces[0])
        fas.append(pieces[1])
    # All done!
    input_file.close()

    return freqs, fas

def read_fas_eas_file(fas_input_file):
    """
    Reads the fas_input_file, returning
    freq, fas_h1, fas_h2, eas, seas
    """
    freqs = []
    fas_h1 = []
    fas_h2 = []
    eas = []
    seas = []

    input_file = open(fas_input_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        if len(pieces) != 5:
            continue
        pieces = [float(piece) for piece in pieces]
        freqs.append(pieces[0])
        fas_h1.append(pieces[1])
        fas_h2.append(pieces[2])
        eas.append(pieces[3])
        seas.append(pieces[4])
    input_file.close()

    return freqs, fas_h1, fas_h2, eas, seas

def plot_fas_comparison(station, sim_file, data_file, label1, label2,
                        outfile, lfreq=None, hfreq=None):
    """
    Plots the FAS comparison between simulated and observed seismograms
    """

    # Set up ticks to match matplotlib 1.x style
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['lines.dashed_pattern'] = [6, 6]
    mpl.rcParams['lines.dashdot_pattern'] = [3, 5, 1, 5]
    mpl.rcParams['lines.dotted_pattern'] = [1, 3]
    mpl.rcParams['lines.scale_dashes'] = False

    # Read data
    (sim_freqs, sim_fas_h1,
     sim_fas_h2, sim_eas, sim_s_eas) = read_fas_eas_file(sim_file)
    (data_freqs, data_fas_h1,
     data_fas_h2, data_eas, data_s_eas) = read_fas_eas_file(data_file)

    # Start plot
    pylab.clf()

    # Figure out limits for x and y axis
    min_x = min([min((sim_freqs)), min(data_freqs)])
    max_x = max([max((sim_freqs)), max(data_freqs)])
    min_horiz_y = min([min(sim_fas_h1), min(data_fas_h1),
                       min(sim_fas_h2), min(data_fas_h2)]) / 1.1
    max_horiz_y = 1.1 * max([max(sim_fas_h1), max(data_fas_h1),
                             max(sim_fas_h2), max(data_fas_h2)])
    min_vert_y = min([min(sim_fas_h1), min(data_fas_h1),
                      min(sim_fas_h2), min(data_fas_h2)]) / 1.1
    max_vert_y = 1.1 * max([max(sim_fas_h1), max(data_fas_h1),
                            max(sim_fas_h2), max(data_fas_h2)])

    pylab.suptitle('FAS for station %s' % (station), size=14)

    pylab.subplots_adjust(top=0.85)
    pylab.subplots_adjust(bottom=0.15)
    pylab.subplots_adjust(left=0.075)
    pylab.subplots_adjust(right=0.975)
    pylab.subplots_adjust(hspace=0.3)
    pylab.subplots_adjust(wspace=0.3)

    # First plot
    ax1 = pylab.subplot(131)
    ax1.set_title('%s, FAS H1' % (label1), fontsize='small')
    pylab.plot(data_freqs, data_fas_h1, label=str(label1),
               linewidth=0.5, color='k')
    pylab.plot(sim_freqs, sim_fas_h1, label=str(label2),
               linewidth=0.5, color='r')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel("Fourier Amplitude")
    pylab.axis([0.01, 100, 0.001, 1000])
    pylab.grid(True)
    pylab.grid(which='major', linestyle='-', color='lightgray')
    #pylab.grid(b=True, which='minor', linewidth=0.5, color='gray')
    if lfreq is not None:
        pylab.vlines(lfreq, 0.001, 1000,
                     color='violet', linestyles='--')
    if hfreq is not None:
        pylab.vlines(hfreq, 0.001, 1000,
                     color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))

    # Second plot
    ax2 = pylab.subplot(132)
    ax2.set_title('%s, FAS H2' % (label1), fontsize='small')
    pylab.plot(data_freqs, data_fas_h2, label=str(label1),
               linewidth=0.5, color='k')
    pylab.plot(sim_freqs, sim_fas_h2, label=str(label2),
               linewidth=0.5, color='r')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel("Fourier Amplitude")
    pylab.axis([0.01, 100, 0.001, 1000])
    pylab.grid(True)
    pylab.grid(which='major', linestyle='-', color='lightgray')
    if lfreq is not None:
        pylab.vlines(lfreq, 0.001, 1000,
                     color='violet', linestyles='--')
    if hfreq is not None:
        pylab.vlines(hfreq, 0.001, 1000,
                     color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))

    # Third plot
    ax3 = pylab.subplot(133)
    ax3.set_title('%s, Smoothed EAS' % (label1), fontsize='small')
    pylab.plot(data_freqs, data_s_eas, label=str(label1),
               linewidth=0.5, color='k')
    pylab.plot(sim_freqs, sim_s_eas, label=str(label2),
               linewidth=0.5, color='r')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel("Fourier Amplitude")
    pylab.axis([0.01, 100, 0.001, 1000])
    pylab.grid(True)
    pylab.grid(which='major', linestyle='-', color='lightgray')
    if lfreq is not None:
        pylab.vlines(lfreq, 0.001, 1000,
                     color='violet', linestyles='--')
    if hfreq is not None:
        pylab.vlines(hfreq, 0.001, 1000,
                     color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))

    pylab.gcf().set_size_inches(10, 4)
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

def plot_fas(freqs, ns_data, ew_data, eas_smoothed_data,
             fas_plot, station, plot_title=None):
    """
    Create a plot of both FAS components
    """

    # Set plot dims
    pylab.gcf().set_size_inches(11, 8.5)
    pylab.gcf().clf()

    # Adjust title y-position
    if plot_title is None:
        plot_title = "Station : %s" % (station)
    t = pylab.title(plot_title, size=12)

    pylab.plot(freqs, ns_data, 'b', lw=0.75, label="NS")
    pylab.plot(freqs, ew_data, 'r', lw=0.75, label="EW")
    pylab.plot(freqs, eas_smoothed_data, 'k', lw=1.25, label="Smoothed EAS")
    pylab.legend(loc='upper right')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.ylabel('Fourier Amplitude (cm/s)')
    pylab.xlabel('Frequency (Hz)')
    pylab.axis([0.01, 100, 0.001, 1000])
    pylab.grid(True)
    pylab.grid(which='major', linestyle='-', color='lightgray')
    pylab.grid(which='minor', linewidth=0.5, color='gray')

    # Save plot
    pylab.savefig(fas_plot, format="png",
                  transparent=False, dpi=plot_config.dpi)
    pylab.close()

def ko98_smoothing(freqs, data, delta_freq, bexp):
    """
    # ** smoothing of a function y (equally-spaced, dx) with the "Konno-Ohmachi"
    # ** function sin (alog10(f/fc)^exp) / alog10(f/fc)^exp) ^^4
    # ** where fc is the frequency around which the smoothing is performed
    # ** exp determines the exponent 10^(1/exp) is the half-width of the peak
    # ** cf Konno & Ohmachi, 1998, BSSA 88-1, pp. 228-241
    """

    nx = len(freqs)
    data_smooth = np.zeros(nx)
    fratio = np.power(10., (2.5 / bexp))
    data_smooth[0] = data[0]

    for index in range(1, nx):
        freq = freqs[index]
        # Added check to avoid division by zero later and NaNs in the output file
        if freq == 0.0:
            data_smooth[index] = data[index]
            continue
        fc1 = freq / fratio
        fc2 = freq * fratio
        index1 = int(fc1 / delta_freq)
        index2 = int((fc2 / delta_freq) + 1)
        if index1 <= 1:
            index1 = 0
        if index2 >= nx:
            index2 = nx
        a1 = 0.0
        a2 = 0.0
        for j in range(index1, index2):
            if j != index:
                # Extra check to avoid NaNs in output file
                if freqs[j] == 0.0:
                    data_smooth[index] = data[index]
                    break
                c1 = bexp * np.log10(freqs[j] / freq)
                c1 = np.power(np.sin(c1) / c1, 4.0)
                a2 = a2 + c1
                a1 = a1 + c1 * data[j]
            else:
                a2 = a2 + 1.0
                a1 = a1 + data[index]
            data_smooth[index] = a1 / a2

    return data_smooth

def calculate_smoothed_eas(ns_file, ew_file, output_file=None):
    """
    Calculates the smoothed EAS at the same frequencies as specified in
    the input files
    """
    b_param = 188.5 # cm/s

    # Read data
    freqs, ns_data = read_fas_file(ns_file)
    _, ew_data = read_fas_file(ew_file)
    eas_data = []

    # Calculate EAS
    for ns_comp, ew_comp in zip(ns_data, ew_data):
        eas_data.append(np.sqrt(0.5*(pow(ns_comp, 2) + pow(ew_comp, 2))))

    # Calculate Smoothed EAS
    smoothed_eas = ko98_smoothing(freqs, eas_data,
                                  freqs[1]-freqs[0],
                                  b_param)

    # Write data file if output_file is provided
    if output_file is not None:
        out_file = open(output_file, 'w')
        out_file.write("# Freq(Hz)\t FAS H1 (cm/s)\t FAS H2 (cm/s)\t "
                       "EAS (cm/s)\t Smoothed EAS, b=%f (cm/s)\n" %
                       (b_param))
        for freq, fas_h1, fas_h2, eas, s_eas in zip(freqs, ns_data,
                                                    ew_data, eas_data,
                                                    smoothed_eas):
            out_file.write("%2.7E\t%2.7E\t%2.7E\t%2.7E\t%2.7E\n" %
                           (freq, fas_h1, fas_h2, eas, s_eas))
        out_file.close()

    # All done!
    return freqs, ns_data, ew_data, eas_data, smoothed_eas

class FAS(object):
    """
    Implement FAS analisys for the Broadband Platform
    """
    def __init__(self, i_r_stations, i_comparison_label=None, sim_id=0):
        """
        Initializes class variables
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.comp_label = i_comparison_label

    def compute_fas(self, a_tmpdir, a_outdir_fas, acc_file,
                    station_name, output_prefix=""):
        """
        Computes FAS for a station
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        asc2smc_control_file = "asc2smc.ctl"
        smc2fs2_control_file = "smc2fs2.ctl"
        header_lines = bband_utils.count_header_lines(os.path.join(a_tmpdir,
                                                                   acc_file))
        # Work on both NS and EW components
        for comp, data_column in zip(["NS", "EW"], [2, 3]):
            # First we convert from BBP to SMC format
            create_boore_asc2smc(os.path.join(a_tmpdir,
                                              asc2smc_control_file),
                                 acc_file, data_column, header_lines,
                                 ".smc8.%s" % (comp))
            cmd = ("%s << END >> %s 2>&1\n" %
                   (os.path.join(install.A_USGS_BIN_DIR, "asc2smc"),
                    self.log) +
                   "%s\n" % (asc2smc_control_file) +
                   "END\n")
            bband_utils.runprog(cmd, False, abort_on_error=True)
            # Then, we run the smc2fs2 FAS tool
            smc_file = "%s.smc8.%s" % (acc_file, comp)
            create_boore_smc2fs2(os.path.join(a_tmpdir,
                                              smc2fs2_control_file),
                                 smc_file, ".no_smooth.fs.col")
            cmd = ("%s >> %s 2>&1\n" %
                   (os.path.join(install.A_USGS_BIN_DIR, "smc2fs2"),
                    self.log))
            bband_utils.runprog(cmd, False, abort_on_error=True)

        # Calculate EAS and smoothed EAS
        ns_file = os.path.join(a_tmpdir,
                               "%s.smc8.NS.no_smooth.fs.col" % (acc_file))
        ew_file = os.path.join(a_tmpdir,
                               "%s.smc8.EW.no_smooth.fs.col" % (acc_file))
        output_file = os.path.join(a_outdir_fas,
                                   "%s%s.smc8.smooth.fs.col" %
                                   (output_prefix, station_name))
        (freqs, ns_fas,
         ew_fas, eas, smoothed_eas) = calculate_smoothed_eas(ns_file,
                                                             ew_file,
                                                             output_file)

        # Create plot
        fas_plot = os.path.join(a_outdir_fas,
                                "%s%s.fas.png" % (output_prefix,
                                                  station_name))
        if output_prefix.startswith("obs"):
            plot_title = "Station: %s - Observed Data" % (station_name)
        else:
            plot_title = None

        plot_fas(freqs, ns_fas, ew_fas, smoothed_eas,
                 fas_plot, station_name, plot_title=plot_title)

    def resample_bbp_file(self, input_bbp_file, output_bbp_file,
                          site, old_dt, new_dt, a_tmpdir):
        """
        Resamples input_bbp_file to new_dt, writing the result to
        output_bbp_file
        """
        sim_id = self.sim_id
        install = install_cfg.InstallCfg.getInstance()

        print("--> Resampling   %.4f --> %.4f - %s... " %
              (old_dt, new_dt, os.path.basename(input_bbp_file)))
        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "%d.%s-fas.acc.000" % (sim_id, site))
        ewfile = os.path.join(a_tmpdir,
                              "%d.%s-fas.acc.090" % (sim_id, site))
        udfile = os.path.join(a_tmpdir,
                              "%d.%s-fas.acc.ver" % (sim_id, site))
        bbpfile = input_bbp_file

        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        progstring = ("%s " %
                      (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                      "nsfile=%s ewfile=%s udfile=%s " %
                      (nsfile, ewfile, udfile) +
                      "wcc2bbp=0 < %s >> %s 2>&1" %
                      (bbpfile, self.log))
        bband_utils.runprog(progstring, abort_on_error=True,
                            print_cmd=False)

        for component in COMPS:
            # Set up paths for each component file
            infile = os.path.join(a_tmpdir,
                                  "%d.%s-fas.acc.%s" %
                                  (sim_id, site, component))

            outfile = os.path.join(a_tmpdir, "%d.%s-fas-resamp.%s" %
                                   (sim_id, site, component))
            bband_utils.check_path_lengths([infile, outfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s newdt=%f " %
                          (os.path.join(install.A_GP_BIN_DIR,
                                        "wcc_resamp_arbdt"), new_dt) +
                          "infile=%s outfile=%s >> %s 2>&1" %
                          (infile, outfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

        # Put the acc file back together
        nsfile = os.path.join(a_tmpdir, "%d.%s-fas-resamp.000" %
                              (sim_id, site))
        ewfile = os.path.join(a_tmpdir, "%d.%s-fas-resamp.090" %
                              (sim_id, site))
        udfile = os.path.join(a_tmpdir, "%d.%s-fas-resamp.ver" %
                              (sim_id, site))
        bbpfile = output_bbp_file
        cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
               "nsfile=%s ewfile=%s udfile=%s " %
               (nsfile, ewfile, udfile) +
               "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
               (bbpfile, self.log))
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

    def run(self):
        """
        Run FAS analysis codes
        """
        print("FAS Calculation".center(80, '-'))

        if self.comp_label is None:
            self.run_scenario()
        else:
            self.run_validation()

        print("FAS Calculation Completed".center(80, '-'))

    def run_scenario(self):
        """
        Run FAS analysis codes for scenario simulations
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.fas_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(a_tmpdir, "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_fas = os.path.join(a_outdir, "FAS")

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_outdir, a_outdir_fas], print_cmd=False)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Save current directory
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir)

        # Initial values
        new_dt = 10000000.0
        new_samples = 0

        # Pre-processing, find lowest dt among all input timeseries
        for site in site_list:
            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)

            # Read DT from simulated file
            input_syn_acc_file = os.path.join(a_outdir, acc_file)
            syn_dt = read_bbp_dt(input_syn_acc_file)
            new_dt = min(new_dt, syn_dt)

        print("==> Target DT set to %.4f" % (new_dt))

        syn_rs_array = []

        print("==> Pre-processing step 1: matching dt...")
        for site in site_list:

            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)

            shutil.copy2(os.path.join(a_outdir, acc_file),
                         os.path.join(a_tmpdir, acc_file))

            # Read DT from input file
            input_syn_acc_file = os.path.join(a_tmpdir, acc_file)
            syn_dt = read_bbp_dt(input_syn_acc_file)

            if syn_dt == new_dt:
                input_rs_syn_acc_file = input_syn_acc_file
            else:
                input_rs_syn_acc_file = os.path.join(a_tmpdir,
                                                     "rs-%s" %
                                                     (acc_file))
                self.resample_bbp_file(input_syn_acc_file,
                                       input_rs_syn_acc_file,
                                       site.scode,
                                       syn_dt, new_dt, a_tmpdir)

            # Keep list of files to be used later
            syn_rs_array.append(input_rs_syn_acc_file)

        print("==> Pre-processing step 2: matching samples...")
        for (site,
             input_rs_syn_acc_file) in zip(site_list,
                                           syn_rs_array):
            # Read samples
            syn_samples = read_bbp_samples(input_rs_syn_acc_file)
            new_samples = max(new_samples, syn_samples)

        print("==> Target samples set to %d" % (new_samples))

        for (site,
             input_rs_syn_acc_file) in zip(site_list,
                                           syn_rs_array):

            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)

            print("==> Processing station: %s" % (site.scode))

            # Make sure file has new_samples data points
            syn_samples = read_bbp_samples(input_rs_syn_acc_file)

            if syn_samples == new_samples:
                input_final_syn_acc_file = input_rs_syn_acc_file
            else:
                input_final_syn_acc_file = os.path.join(a_tmpdir,
                                                        "pad-%s" %
                                                        (acc_file))
                extra_points = new_samples - syn_samples
                print("=====> Adding %d extra points to %s" %
                      (extra_points, os.path.basename(input_rs_syn_acc_file)))
                add_extra_points(input_rs_syn_acc_file,
                                 input_final_syn_acc_file,
                                 extra_points)

            # Copy calculated acc file to tmpdata
            print("=====> computing FAS for calculated seismogram... ", end="", flush=True)
            t1 = time.time()
            self.compute_fas(a_tmpdir, a_outdir_fas,
                             os.path.basename(input_final_syn_acc_file),
                             site.scode, output_prefix="%d." % (self.sim_id))
            t2 = time.time()
            print("%10.2f s" % (t2 - t1))

        # All done, restore working directory
        os.chdir(old_cwd)

    def run_validation(self):
        """
        Run FAS analysis codes for validation simulations
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.fas_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(a_tmpdir, "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_fas = os.path.join(a_outdir, "FAS")

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_outdir, a_outdir_fas], print_cmd=False)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Save current directory
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir)

        # Initial values
        new_dt = 10000000.0
        new_samples = 0

        # Pre-processing, find lowest dt among all input timeseries
        for site in site_list:
            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)
            obs_acc_file = "%s.bbp" % (site.scode)

            # Read DT from both files
            input_syn_acc_file = os.path.join(a_outdir, acc_file)
            input_obs_acc_file = os.path.join(a_tmpdir_seis, obs_acc_file)
            syn_dt = read_bbp_dt(input_syn_acc_file)
            obs_dt = read_bbp_dt(input_obs_acc_file)
            new_station_dt = min(syn_dt, obs_dt)
            new_dt = min(new_dt, new_station_dt)

        print("==> Target DT set to %.4f" % (new_dt))

        obs_rs_array = []
        syn_rs_array = []

        print("==> Pre-processing step 1: matching dt...")
        for site in site_list:

            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)
            obs_acc_file = "%s.bbp" % (site.scode)

            shutil.copy2(os.path.join(a_outdir, acc_file),
                         os.path.join(a_tmpdir, acc_file))
            shutil.copy2(os.path.join(a_tmpdir_seis, obs_acc_file),
                         os.path.join(a_tmpdir, obs_acc_file))

            # Read DT from both files
            input_syn_acc_file = os.path.join(a_tmpdir, acc_file)
            input_obs_acc_file = os.path.join(a_tmpdir, obs_acc_file)
            syn_dt = read_bbp_dt(input_syn_acc_file)
            obs_dt = read_bbp_dt(input_obs_acc_file)

            if syn_dt == new_dt:
                input_rs_syn_acc_file = input_syn_acc_file
            else:
                input_rs_syn_acc_file = os.path.join(a_tmpdir,
                                                     "rs-%s" %
                                                     (acc_file))
                self.resample_bbp_file(input_syn_acc_file,
                                       input_rs_syn_acc_file,
                                       site.scode,
                                       syn_dt, new_dt, a_tmpdir)
            if obs_dt == new_dt:
                input_rs_obs_acc_file = input_obs_acc_file
            else:
                input_rs_obs_acc_file = os.path.join(a_tmpdir,
                                                     "rs-%s" %
                                                     (obs_acc_file))
                self.resample_bbp_file(input_obs_acc_file,
                                       input_rs_obs_acc_file,
                                       site.scode,
                                       obs_dt, new_dt, a_tmpdir)

            # Keep list of files to be used later
            syn_rs_array.append(input_rs_syn_acc_file)
            obs_rs_array.append(input_rs_obs_acc_file)

        print("==> Pre-processing step 2: matching samples...")
        for (site,
             input_rs_syn_acc_file,
             input_rs_obs_acc_file) in zip(site_list,
                                           syn_rs_array,
                                           obs_rs_array):
            # Read samples
            syn_samples = read_bbp_samples(input_rs_syn_acc_file)
            obs_samples = read_bbp_samples(input_rs_obs_acc_file)
            new_station_samples = max(syn_samples, obs_samples)
            new_samples = max(new_samples, new_station_samples)

        print("==> Target samples set to %d" % (new_samples))

        for (site,
             input_rs_syn_acc_file,
             input_rs_obs_acc_file) in zip(site_list,
                                           syn_rs_array,
                                           obs_rs_array):

            acc_file = "%d.%s.acc.bbp" % (sim_id, site.scode)
            obs_acc_file = "%s.bbp" % (site.scode)

            print("==> Processing station: %s" % (site.scode))

            # Make sure both files have equal number of data points
            syn_samples = read_bbp_samples(input_rs_syn_acc_file)
            obs_samples = read_bbp_samples(input_rs_obs_acc_file)

            if syn_samples == new_samples:
                input_final_syn_acc_file = input_rs_syn_acc_file
            else:
                input_final_syn_acc_file = os.path.join(a_tmpdir,
                                                        "pad-%s" %
                                                        (acc_file))
                extra_points = new_samples - syn_samples
                print("=====> Adding %d extra points to %s" %
                      (extra_points, os.path.basename(input_rs_syn_acc_file)))
                add_extra_points(input_rs_syn_acc_file,
                                 input_final_syn_acc_file,
                                 extra_points)

            if obs_samples == new_samples:
                input_final_obs_acc_file = input_rs_obs_acc_file
            else:
                extra_points = new_samples - obs_samples
                input_final_obs_acc_file = os.path.join(a_tmpdir,
                                                        "pad-%s" %
                                                        (obs_acc_file))
                print("=====> Adding %d extra points to %s" %
                      (extra_points, os.path.basename(input_rs_obs_acc_file)))
                add_extra_points(input_rs_obs_acc_file,
                                 input_final_obs_acc_file,
                                 extra_points)

            # Copy calculated acc file to tmpdata
            print("=====> computing FAS for calculated seismogram... ", end="", flush=True)
            t1 = time.time()
            self.compute_fas(a_tmpdir, a_outdir_fas,
                             os.path.basename(input_final_syn_acc_file),
                             site.scode, output_prefix="%d." % (self.sim_id))
            t2 = time.time()
            print("%10.2f s" % (t2 - t1))


            print("=====> computing FAS for observed seismogram... ", end="", flush=True)
            t1 = time.time()
            self.compute_fas(a_tmpdir, a_outdir_fas,
                             os.path.basename(input_final_obs_acc_file),
                             site.scode, output_prefix="obs.")
            t2 = time.time()
            print("%10.2f s" % (t2 - t1))

            # Create comparison plot
            max_syn_freq = 1.0 / (syn_dt * 2)
            hfreq = min(max_syn_freq, site.high_freq_corner)
            sim_file = os.path.join(a_outdir_fas,
                                    "%d.%s.smc8.smooth.fs.col" %
                                    (self.sim_id, site.scode))
            data_file = os.path.join(a_outdir_fas,
                                     "obs.%s.smc8.smooth.fs.col" % (site.scode))
            outfile = os.path.join(a_outdir_fas,
                                   "%d.%s.fas.comparison.png" %
                                   (self.sim_id, site.scode))
            plot_fas_comparison(site.scode, sim_file, data_file,
                                self.comp_label, self.sim_id, outfile,
                                lfreq=site.low_freq_corner, hfreq=hfreq)

        # All done, restore working directory
        os.chdir(old_cwd)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: %s station_list comparison_label sim_id" %
              (os.path.basename(sys.argv[0])))
        sys.exit(1)
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = FAS(sys.argv[1], sys.argv[2], sim_id=int(sys.argv[3]))
    ME.run()
    sys.exit(0)
