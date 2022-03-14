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

Plots Goodness of Fit bias +/- sigma
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import pylab
import matplotlib.gridspec as gridspec

# Import plot config file
import bband_utils
import plot_config

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
MIN_Y_AXIS_RATIO = -0.5
MAX_Y_AXIS_RATIO = 0.5
MAX_PERIOD = 10.0
XTICK_LOC = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
XTICK_LABEL = ['0.1', '0.2', '0.5', '1', '2', '5', '10']
XTICK_LOC_0_01 = [0.01, 0.02, 0.05,
                  0.1, 0.2, 0.5,
                  1.0, 2.0, 5.0, 10.0]
XTICK_LABEL_0_01 = ['0.01', '0.02', '0.05',
                    '0.1', '0.2', '0.5',
                    '1', '2', '5', '10']
XTICK_FREQ_LOC = [10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1]
XTICK_FREQ_LABEL = ['10', '5', '2', '1', '0.5', '0.2', '0.1']
XTICK_FREQ_LOC_0_01 = [100.0, 50.0, 20.0, 10.0,
                       5.0, 2.0, 1.0, 0.5, 0.2,
                       0.1, 0.05, 0.02, 0.01]
XTICK_FREQ_LABEL_0_01 = ['100', '50', '20',
                         '10', '5', '2',
                         '1', '0.5', '0.2',
                         '0.1', '0.05', '0.02', '0.01']
XTICK_FAS_LOC_0_01 = [0.01, 0.02, 0.05,
                      0.1, 0.2, 0.5,
                      1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
XTICK_FAS_LABEL_0_01 = ['0.01', '0.02', '0.05',
                        '0.1', '0.2', '0.5',
                        '1', '2', '5', '10',
                        '20', '50', '100']

# Component extentions
COMP_EXT_RD50 = ['psa5n', 'psa5e', 'rotd50']
COMP_EXT_RD100 = ['rotd100', 'rotd50', 'ratio']
COMP_EXT_FAS = ['fash1', 'fash2', 'seas']

# Component titles
COMP_TITLE_RD50 = ['PSA North 5%', 'PSA East 5%', 'RotD50']
COMP_TITLE_RD100 = ['RotD100', 'RotD50', 'RotD100/RotD50']
COMP_TITLE_FAS = ['FAS North', 'FAS East', 'SEAS']

# Component subplot locations
COMP_OFFSET = [312, 313, 311]
#COMP_OFFSET_FAS = [211, 212, 222]
COMP_OFFSET_FAS = [312, 313, 311]

# Colorsets for plots
COLORSETS = {"single": 0,
             "combined": 1}
BIAS_COLORS = ["red", "black"]
BIAS_LH_COLORS = ["cyan", (0.8, 0.8, 1)]
CONF_LH_COLORS = ["yellow", (1, 0.3, 0.9)]

class PlotGoF(object):
    def __init__(self):
        return

    def read_data(self, datafile, min_period, max_period=MAX_PERIOD):
        """
        Read in response spectra data from specified datafile
        """
        data = [[], [],]

        # Read input file
        in_file = open(datafile, "r")
        for line in in_file:
            if line.startswith("#") or line.startswith("%"):
                continue
            tmp = line.split()
            period = float(tmp[0])
            # Extract subset of period values
            if ((period >= min_period) and
                (period <= max_period)):
                data[0].append(float(tmp[0]))
                data[1].append(float(tmp[1]))
        # Close file
        in_file.close()
        # Return data
        return data

    def multi_plot(self, plottitle, gof_fileroot, indir,
                   outdir, legends, num_stations, mode="P"):
        """
        Creates several GOF plots using the files specified in the
        gof_fileroot array. mode selects periods (P) or frequencies (F)
        for the X axis.
        """
        mode = mode.upper()
        if mode != "P" and mode != "F":
            print("Invalid mode, must specify 'P' for periods or 'F' for frequencies!")
            return
        # Pick components and labels
        if mode == "P":
            xtick_loc = XTICK_LOC_0_01
            xtick_label = XTICK_LABEL_0_01
        elif mode == "F":
            xtick_loc = XTICK_FREQ_LOC_0_01
            xtick_label = XTICK_FREQ_LABEL_0_01

        # Initialize data arrays
        freqs = [[] for _ in range(len(gof_fileroot))]
        periods = [[] for _ in range(len(gof_fileroot))]
        bias = [[] for _ in range(len(gof_fileroot))]
        m90 = [[] for _ in range(len(gof_fileroot))]
        p90 = [[] for _ in range(len(gof_fileroot))]
        sigma = [[] for _ in range(len(gof_fileroot))]
        sigma0 = [[] for _ in range(len(gof_fileroot))]

        bias_l = [[] for _ in range(len(gof_fileroot))]
        bias_h = [[] for _ in range(len(gof_fileroot))]
        conf_l = [[] for _ in range(len(gof_fileroot))]
        conf_h = [[] for _ in range(len(gof_fileroot))]

        # Read data from all input files
        for compnum in range(0, len(gof_fileroot)):
            filenamebase = os.path.join(indir, gof_fileroot[compnum])
            #print("Reading component files %s.*" % (filenamebase))
            periods[compnum], bias[compnum] = self.read_data("%s.bias" %
                                                             (filenamebase),
                                                             0.01)
            periods[compnum], m90[compnum] = self.read_data("%s.m90" %
                                                            (filenamebase),
                                                            0.01)
            periods[compnum], p90[compnum] = self.read_data("%s.p90" %
                                                            (filenamebase),
                                                            0.01)
            periods[compnum], sigma[compnum] = self.read_data("%s.sigma" %
                                                              (filenamebase),
                                                              0.01)
            periods[compnum], sigma0[compnum] = self.read_data("%s.sigma0" %
                                                               (filenamebase),
                                                               0.01)

            # Compute bias and conf interval lower/upper bounds
            for i in range(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

            if mode == "F":
                # Add extra point for T=100s
                for i in range(0, len(bias[compnum])):
                    bias[compnum].append(0.0)
                    bias_l[compnum].append(0.0)
                    bias_h[compnum].append(0.0)
                    conf_l[compnum].append(0.0)
                    conf_h[compnum].append(0.0)
                    bias[compnum].append(0.0)
                    bias_l[compnum].append(0.0)
                    bias_h[compnum].append(0.0)
                    conf_l[compnum].append(0.0)
                    conf_h[compnum].append(0.0)
                    periods[compnum].append(10.01)
                    periods[compnum].append(100.0)

        num_periods = len(periods[0])
        for comp in periods:
            if len(comp) != num_periods:
                print("Number of data points unequal across components")
                return

        # Calculate frequencies
        for freq, period in zip(freqs, periods):
            for item in period:
                freq.append(1.0 / item)

        # Construct baseline
        baseline = [0.0 for _ in periods[0]]

        # Find max, min values
        if mode == "P":
            min_x = min([min(comp) for comp in periods])
            max_x = max([max(comp) for comp in periods])
        elif mode == "F":
            min_x = min(freqs[0])
            max_x = max(freqs[0])
            for comp in freqs:
                min_x = min(min_x, min(comp))
                max_x = max(max_x, max(comp))
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Set up ticks to match matplotlib 1.x style
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True

        # Start plots
        num_plots = len(gof_fileroot)
        # Make 1 column of "num_plots" rows
        fig, axs = pylab.plt.subplots(num_plots, 1)
        # Set plot dims
        fig.set_size_inches(6, 10)
        fig.subplots_adjust(left=0.1)
        fig.subplots_adjust(right=0.97)
        fig.subplots_adjust(top=0.92)
        fig.subplots_adjust(bottom=0.07)
        fig.subplots_adjust(hspace=0.4)
        fig.subplots_adjust(wspace=0.5)

        # Add subplots in a list
        subfigs = []
        for idx in range(0, num_plots):
            subfigs.append(axs[idx])

        # Now walk through each subfig
        for (subfig, subplot_title,
             cur_period, cur_freq,
             cur_bias, cur_bias_h, cur_bias_l,
             cur_conf_h, cur_conf_l) in zip(subfigs,
                                            legends,
                                            periods,
                                            freqs,
                                            bias,
                                            bias_h,
                                            bias_l,
                                            conf_h,
                                            conf_l):
            subfig.set_xlim(min_x, max_x)
            subfig.set_ylim(min_y, max_y)
            if mode == "P":
                x_comp = cur_period
            elif mode == "F":
                x_comp = cur_freq
            subfig.set_title("%s" % subplot_title, size=10)
            subfig.plot(x_comp, cur_bias, color='red',
                        label='_nolegend_', linewidth=1.0)
            subfig.fill_between(x_comp, cur_bias_h,
                                cur_bias_l, color='cyan',
                                label='_nolegend_')
            subfig.fill_between(x_comp, cur_conf_h,
                                cur_conf_l, color='yellow',
                                label='_nolegend_')
            subfig.plot(x_comp, baseline, color='grey',
                        label='_nolegend_', linewidth=1.0)
            if mode == "P":
                subfig.set_xlabel("Period (sec)", size=8)
            elif mode == "F":
                subfig.set_xlabel("Frequency (Hz)", size=8)
            subfig.set_ylabel("ln (data/model)", size=8)
            subfig.set_xscale('log')
            # Old way to do it
            # subfig.set_xticks(xtick_loc, xtick_label)
            subfig.set_xticks(xtick_loc)
            subfig.set_xticklabels(xtick_label)
            subfig.tick_params(labelsize=8)
            subfig.minorticks_on()

        fig.suptitle('%s\nNumber of stations: %d' % (plottitle, num_stations),
                     size=12)
        # Figure out output filename
        outfile = gof_fileroot[0]
        outfile = outfile[:outfile.rfind("-")]
        if mode == "P":
            outfile = os.path.join(outdir, "gof-%s.png" % (outfile))
        elif mode == "F":
            outfile = os.path.join(outdir, "gof-%s-freq.png" % (outfile))
        print("==> Created GoF plot: %s" % (outfile))
        fig.savefig(outfile, format="png",
                    transparent=False, dpi=plot_config.dpi)
        pylab.close()

    def plot_single_component_gof(self, plottitle, gof_fileroot, indir, outdir,
                                  cutoff=0, min_period=0.01, colorset=None):
        """
        Creates a single component GOF plot (e.g. RotD50)
        """
        # Pick components and labels
        if min_period == 0.01:
            xtick_loc = XTICK_LOC_0_01
            xtick_label = XTICK_LABEL_0_01
        elif min_period == 0.1:
            xtick_loc = XTICK_LOC
            xtick_label = XTICK_LABEL
        else:
            raise bband_utils.ParameterError("invalid min_period: %f" %
                                             (min_period))

        comp_ext = COMP_EXT_RD50

        # Select colorset
        if colorset is None:
            colorset_idx = 0
        else:
            colorset_idx = COLORSETS[colorset]

        # Set up ticks to match matplotlib 1.x style
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True

        # Set plot dims
        pylab.gcf().set_size_inches(6, 2.8)
        pylab.gcf().clf()

        pylab.subplots_adjust(left=0.10)
        pylab.subplots_adjust(right=0.97)
        pylab.subplots_adjust(top=0.77)
        pylab.subplots_adjust(bottom=0.13)
        pylab.subplots_adjust(hspace=0.35)
        pylab.subplots_adjust(wspace=0.5)

        period = [[], [], []]
        bias = [[], [], []]
        m90 = [[], [], []]
        p90 = [[], [], []]
        sigma = [[], [], []]
        sigma0 = [[], [], []]

        bias_l = [[], [], []]
        bias_h = [[], [], []]
        conf_l = [[], [], []]
        conf_h = [[], [], []]

        for compnum in range(0, len(comp_ext)):
            comp = comp_ext[compnum]
            filenamebase = os.path.join(indir, "%s-%s" % (gof_fileroot, comp))
            #print("Reading component files %s.*" % (filenamebase))
            period[compnum], bias[compnum] = self.read_data("%s.bias" %
                                                            (filenamebase),
                                                            min_period)
            period[compnum], m90[compnum] = self.read_data("%s.m90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], p90[compnum] = self.read_data("%s.p90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], sigma[compnum] = self.read_data("%s.sigma" %
                                                             (filenamebase),
                                                             min_period)
            period[compnum], sigma0[compnum] = self.read_data("%s.sigma0" %
                                                              (filenamebase),
                                                              min_period)

            # Compute bias and conf interval lower/upper bounds
            for i in range(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

        # Make sure all components have same number of data points
        npts = [len(component) for component in period]
        if npts[1:] != npts[:-1]:
            print("Number of data points unequal across components")
            return

        # Construct baseline
        baseline = []
        for _ in range(0, len(period[0])):
            baseline.append(0.0)

        # Find max, min values
        min_x = min(period[0])
        max_x = max(period[0])
        for comp in period:
            min_x = min(min_x, min(comp))
            max_x = max(max_x, max(comp))
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Draw each component
        compnum = 2
        comp = 'rotd50'

        pylab.subplot(111)
        pylab.title('RotD50', size='small')
        pylab.plot(period[compnum], bias[compnum],
                   color=BIAS_COLORS[colorset_idx],
                   label='_nolegend_', linewidth=1.0)
        pylab.fill_between(period[compnum], bias_h[compnum],
                           bias_l[compnum],
                           color=BIAS_LH_COLORS[colorset_idx],
                           label='_nolegend_')
        pylab.fill_between(period[compnum], conf_h[compnum],
                           conf_l[compnum],
                           color=CONF_LH_COLORS[colorset_idx],
                           label='_nolegend_')
        pylab.plot(period[compnum], baseline, color='grey',
                   label='_nolegend_', linewidth=1.0)
        pylab.xlim(min_x, max_x)
        pylab.ylim(min_y, max_y)
        pylab.xlabel("Period (sec)", size=8)
        pylab.ylabel("ln (data/model)", size=8)
        pylab.xscale('log')
        pylab.xticks(xtick_loc, xtick_label)
        pylab.tick_params(labelsize=8)

        if cutoff == 0:
            pylab.suptitle('%s' % (plottitle), size=11)
        else:
            pylab.suptitle('%s\nR < %d km' % (plottitle, cutoff), size=11)
        outfile = os.path.join(outdir, "gof-%s.png" % (gof_fileroot))
        print("==> Created GoF plot: %s" % (outfile))
        pylab.savefig(outfile, format="png",
                      transparent=False, dpi=plot_config.dpi)
        pylab.close()

    def plot_single_component_freq_gof(self, plottitle, gof_fileroot, indir, outdir,
                                       cutoff=0, min_period=0.01, colorset=None):
        """
        Creates a single component GOF plot using frequencies instead of periods
        """
        # Pick components and labels
        if min_period == 0.01:
            xtick_loc = XTICK_FREQ_LOC_0_01
            xtick_label = XTICK_FREQ_LABEL_0_01
        elif min_period == 0.1:
            xtick_loc = XTICK_FREQ_LOC
            xtick_label = XTICK_FREQ_LABEL
        else:
            raise bband_utils.ParameterError("invalid min_period: %f" %
                                             (min_period))

        comp_ext = COMP_EXT_RD50

        # Select colorset
        if colorset is None:
            colorset_idx = 0
        else:
            colorset_idx = COLORSETS[colorset]

        # Set up ticks to match matplotlib 1.x style
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True

        # Set plot dims
        pylab.gcf().set_size_inches(6, 2.8)
        pylab.gcf().clf()

        pylab.subplots_adjust(left=0.10)
        pylab.subplots_adjust(right=0.97)
        pylab.subplots_adjust(top=0.77)
        pylab.subplots_adjust(bottom=0.13)
        pylab.subplots_adjust(hspace=0.35)
        pylab.subplots_adjust(wspace=0.5)

        freq = [[], [], []]
        period = [[], [], []]
        bias = [[], [], []]
        m90 = [[], [], []]
        p90 = [[], [], []]
        sigma = [[], [], []]
        sigma0 = [[], [], []]

        bias_l = [[], [], []]
        bias_h = [[], [], []]
        conf_l = [[], [], []]
        conf_h = [[], [], []]

        for compnum in range(0, len(comp_ext)):
            comp = comp_ext[compnum]
            filenamebase = os.path.join(indir, "%s-%s" % (gof_fileroot, comp))
            #print("Reading component files %s.*" % (filenamebase))
            period[compnum], bias[compnum] = self.read_data("%s.bias" %
                                                            (filenamebase),
                                                            min_period)
            period[compnum], m90[compnum] = self.read_data("%s.m90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], p90[compnum] = self.read_data("%s.p90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], sigma[compnum] = self.read_data("%s.sigma" %
                                                             (filenamebase),
                                                             min_period)
            period[compnum], sigma0[compnum] = self.read_data("%s.sigma0" %
                                                              (filenamebase),
                                                              min_period)

            # Compute bias and conf interval lower/upper bounds
            for i in range(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

        # Add extra point for T=100s
        for compnum in range(0, len(comp_ext)):
            bias[compnum].append(0.0)
            bias_l[compnum].append(0.0)
            bias_h[compnum].append(0.0)
            conf_l[compnum].append(0.0)
            conf_h[compnum].append(0.0)
            bias[compnum].append(0.0)
            bias_l[compnum].append(0.0)
            bias_h[compnum].append(0.0)
            conf_l[compnum].append(0.0)
            conf_h[compnum].append(0.0)
            period[compnum].append(10.01)
            period[compnum].append(100.0)

        # Calculate frequencies
        for freqs, periods in zip(freq, period):
            for item in periods:
                freqs.append(1.0 / item)

        # Make sure all components have same number of data points
        npts = [len(component) for component in freq]
        if npts[1:] != npts[:-1]:
            print("Number of data points unequal across components")
            return

        # Construct baseline
        baseline = []
        for _ in range(0, len(freq[0])):
            baseline.append(0.0)

        # Find max, min values
        min_x = min(freq[0])
        max_x = max(freq[0])
        for comp in freq:
            min_x = min(min_x, min(comp))
            max_x = max(max_x, max(comp))
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Draw each component
        compnum = 2
        comp = 'rotd50'

        pylab.subplot(111)
        pylab.title('RotD50 PSA', size='small')
        pylab.plot(freq[compnum], bias[compnum],
                   color=BIAS_COLORS[colorset_idx],
                   label='_nolegend_', linewidth=1.0)
        pylab.fill_between(freq[compnum], bias_h[compnum],
                           bias_l[compnum],
                           color=BIAS_LH_COLORS[colorset_idx],
                           label='_nolegend_')
        pylab.fill_between(freq[compnum], conf_h[compnum],
                           conf_l[compnum],
                           color=CONF_LH_COLORS[colorset_idx],
                           label='_nolegend_')
        pylab.plot(freq[compnum], baseline, color='grey',
                   label='_nolegend_', linewidth=1.0)
        pylab.xlim(min_x, max_x)
        pylab.ylim(min_y, max_y)
        pylab.xlabel("Frequency (Hz)", size=8)
        pylab.ylabel("ln (data/model)", size=8)
        pylab.xscale('log')
        pylab.xticks(xtick_loc, xtick_label)
        pylab.tick_params(labelsize=8)

        if cutoff == 0:
            pylab.suptitle('%s' % (plottitle), size=11)
        else:
            pylab.suptitle('%s\nR < %d km' % (plottitle, cutoff), size=11)
        outfile = os.path.join(outdir, "gof-%s-freq.png" % (gof_fileroot))
        print("==> Created GoF plot: %s" % (outfile))
        pylab.savefig(outfile, format="png",
                      transparent=False, dpi=plot_config.dpi)
        pylab.close()

    def plot_three_component_gof(self, plottitle, gof_fileroot, indir, outdir,
                                 cutoff=0, min_period=0.01, mode=None, colorset=None):
        """
        Creates a GOF plot with three subplots (e.g. RotD50/PSA5n/PSA5e)
        """
        # Pick components and labels
        if min_period == 0.01:
            xtick_loc = XTICK_LOC_0_01
            xtick_label = XTICK_LABEL_0_01
        elif min_period == 0.1:
            xtick_loc = XTICK_LOC
            xtick_label = XTICK_LABEL
        else:
            raise bband_utils.ParameterError("invalid min_period: %f" %
                                             (min_period))

        if mode == "rd50":
            comp_ext = COMP_EXT_RD50
            comp_title = COMP_TITLE_RD50
        elif mode == "rd100":
            comp_ext = COMP_EXT_RD100
            comp_title = COMP_TITLE_RD100
        else:
            raise bband_utils.ParameterError("plot mode %s unsupported" %
                                             (mode))

        # Select colorset
        if colorset is None:
            colorset_idx = 0
        else:
            colorset_idx = COLORSETS[colorset]

        # Set up ticks to match matplotlib 1.x style
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True

        # Set plot dims
        pylab.gcf().set_size_inches(6, 8)
        pylab.gcf().clf()

        pylab.subplots_adjust(left=0.125)
        pylab.subplots_adjust(right=0.95)
        pylab.subplots_adjust(top=0.9)
        pylab.subplots_adjust(bottom=0.1)
        pylab.subplots_adjust(hspace=0.35)
        pylab.subplots_adjust(wspace=0.5)

        period = [[], [], []]
        bias = [[], [], []]
        m90 = [[], [], []]
        p90 = [[], [], []]
        sigma = [[], [], []]
        sigma0 = [[], [], []]

        bias_l = [[], [], []]
        bias_h = [[], [], []]
        conf_l = [[], [], []]
        conf_h = [[], [], []]

        for compnum in range(0, len(comp_ext)):
            comp = comp_ext[compnum]
            filenamebase = os.path.join(indir, "%s-%s" % (gof_fileroot, comp))
            #print("Reading component files %s.*" % (filenamebase))
            period[compnum], bias[compnum] = self.read_data("%s.bias" %
                                                            (filenamebase),
                                                            min_period)
            period[compnum], m90[compnum] = self.read_data("%s.m90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], p90[compnum] = self.read_data("%s.p90" %
                                                           (filenamebase),
                                                           min_period)
            period[compnum], sigma[compnum] = self.read_data("%s.sigma" %
                                                             (filenamebase),
                                                             min_period)
            period[compnum], sigma0[compnum] = self.read_data("%s.sigma0" %
                                                              (filenamebase),
                                                              min_period)

            # Compute bias and conf interval lower/upper bounds
            for i in range(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

        # Make sure all components have same number of data points
        npts = [len(component) for component in period]
        if npts[1:] != npts[:-1]:
            print("Number of data points unequal across components")
            return

        # Construct baseline
        baseline = []
        for _ in range(0, len(period[0])):
            baseline.append(0.0)

        # Find max, min values
        min_x = min(period[0])
        max_x = max(period[0])
        for comp in period:
            min_x = min(min_x, min(comp))
            max_x = max(max_x, max(comp))
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Draw each component
        for compnum in range(0, 3):
            comp = comp_ext[compnum]
            offset = COMP_OFFSET[compnum]

            pylab.subplot(offset)
            pylab.title(comp_title[compnum], size='small')
            pylab.plot(period[compnum], bias[compnum],
                       color=BIAS_COLORS[colorset_idx],
                       label='_nolegend_', linewidth=1.0)
            pylab.fill_between(period[compnum], bias_h[compnum],
                               bias_l[compnum],
                               color=BIAS_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.fill_between(period[compnum], conf_h[compnum],
                               conf_l[compnum],
                               color=CONF_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.plot(period[compnum], baseline, color='grey',
                       label='_nolegend_', linewidth=1.0)
            pylab.xlim(min_x, max_x)
            if comp == 'ratio':
                pylab.ylim(MIN_Y_AXIS_RATIO, MAX_Y_AXIS_RATIO)
            else:
                pylab.ylim(min_y, max_y)
            pylab.xlabel("Period (sec)", size=8)
            pylab.ylabel("ln (data/model)", size=8)
            pylab.xscale('log')
            pylab.xticks(xtick_loc, xtick_label)
            pylab.tick_params(labelsize=8)

        if cutoff == 0:
            pylab.suptitle('%s' % (plottitle), size=11)
        else:
            pylab.suptitle('%s\nR < %d km' % (plottitle, cutoff), size=11)
        outfile = os.path.join(outdir, "gof-%s.png" % (gof_fileroot))
        print("==> Created GoF plot: %s" % (outfile))
        pylab.savefig(outfile, format="png",
                      transparent=False, dpi=plot_config.dpi)
        pylab.close()

    def plot(self, plottitle, gof_fileroot, indir, outdir,
             cutoff=0, min_period=0.01, mode=None, colorset=None):
        """
        Creates the GOF plot
        """
        if mode == "rd50-single-freq":
            self.plot_single_component_freq_gof(plottitle, gof_fileroot,
                                                indir, outdir,
                                                cutoff=cutoff,
                                                min_period=min_period,
                                                colorset=colorset)
        elif mode == "rd50-single":
            self.plot_single_component_gof(plottitle, gof_fileroot,
                                           indir, outdir,
                                           cutoff=cutoff,
                                           min_period=min_period,
                                           colorset=colorset)
        elif mode == "rd50" or mode == "rd100":
            self.plot_three_component_gof(plottitle, gof_fileroot,
                                          indir, outdir,
                                          cutoff=cutoff,
                                          min_period=min_period,
                                          mode=mode, colorset=colorset)
        else:
            raise bband_utils.ParameterError("plot mode %s unsupported" %
                                             (mode))

    def plot_fas_gof(self, plottitle, gof_fileroot, indir,
                     outdir, cutoff=0, colorset=None):
        """
        Creates a FAS GOF plot with three subplots (SEAS, FAS_H1, FAS_H2)
        """
        # Pick components, labels, and colorset
        xtick_loc = XTICK_FAS_LOC_0_01
        xtick_label = XTICK_FAS_LABEL_0_01
        comp_ext = COMP_EXT_FAS
        comp_title = COMP_TITLE_FAS
        min_period = 0.01
        max_period = 100.0

        # Select colorset
        if colorset is None:
            colorset_idx = 0
        else:
            colorset_idx = COLORSETS[colorset]

        # Set up ticks to match matplotlib 1.x style
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.top'] = True
        mpl.rcParams['ytick.right'] = True

        # Set plot dims
        pylab.gcf().set_size_inches(6, 9)
        pylab.gcf().clf()

#        # Create 2x2 sub plots
#        gs = gridspec.GridSpec(2, 2)
#        COMP_OFFSET_FAS = [gs[1, 0], gs[1, 1], gs[0, :]]

        pylab.subplots_adjust(left=0.10)
        pylab.subplots_adjust(right=0.97)
        pylab.subplots_adjust(top=0.9)
        pylab.subplots_adjust(bottom=0.1)
        pylab.subplots_adjust(hspace=0.35)
        pylab.subplots_adjust(wspace=0.5)

        period = [[], [], []]
        bias = [[], [], []]
        m90 = [[], [], []]
        p90 = [[], [], []]
        sigma = [[], [], []]
        sigma0 = [[], [], []]

        bias_l = [[], [], []]
        bias_h = [[], [], []]
        conf_l = [[], [], []]
        conf_h = [[], [], []]

        for compnum in range(0, len(comp_ext)):
            comp = comp_ext[compnum]
            filenamebase = os.path.join(indir, "%s-%s" % (gof_fileroot, comp))
            # print("Reading component files %s.*" % (filenamebase))
            period[compnum], bias[compnum] = self.read_data("%s.bias" %
                                                            (filenamebase),
                                                            min_period,
                                                            max_period)
            period[compnum], m90[compnum] = self.read_data("%s.m90" %
                                                           (filenamebase),
                                                           min_period,
                                                           max_period)
            period[compnum], p90[compnum] = self.read_data("%s.p90" %
                                                           (filenamebase),
                                                           min_period,
                                                           max_period)
            period[compnum], sigma[compnum] = self.read_data("%s.sigma" %
                                                             (filenamebase),
                                                             min_period,
                                                             max_period)
            period[compnum], sigma0[compnum] = self.read_data("%s.sigma0" %
                                                              (filenamebase),
                                                              min_period,
                                                              max_period)

            # Compute bias and conf interval lower/upper bounds
            for i in range(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

        # Make sure all components have same number of data points
        npts = [len(component) for component in period]
        #print(npts)
        if npts[1:] != npts[:-1]:
            print("Number of data points unequal across components")
            return

        # Construct baseline
        baseline = []
        for _ in range(0, len(period[0])):
            baseline.append(0.0)

        # Find max, min values
        min_x = min(period[0])
        max_x = max(period[0])
        for comp in period:
            min_x = min(min_x, min(comp))
            max_x = max(max_x, max(comp))
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Draw each component
        for compnum in range(0, 3):
            comp = comp_ext[compnum]
            offset = COMP_OFFSET_FAS[compnum]

            pylab.subplot(offset)
            pylab.title(comp_title[compnum], size='small')
            pylab.plot(period[compnum], bias[compnum],
                       color=BIAS_COLORS[colorset_idx],
                       label='_nolegend_', linewidth=1.0)
            pylab.fill_between(period[compnum], bias_h[compnum],
                               bias_l[compnum],
                               color=BIAS_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.fill_between(period[compnum], conf_h[compnum],
                               conf_l[compnum],
                               color=CONF_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.plot(period[compnum], baseline, color='grey',
                       label='_nolegend_', linewidth=1.0)
            pylab.xlim(min_x, max_x)
            if comp == 'ratio':
                pylab.ylim(MIN_Y_AXIS_RATIO, MAX_Y_AXIS_RATIO)
            else:
                pylab.ylim(min_y, max_y)
            pylab.xlabel("Frequency (Hz)", size=8)
            pylab.ylabel("ln (data/model)", size=8)
            pylab.xscale('log')
            pylab.xticks(xtick_loc, xtick_label)
            pylab.tick_params(labelsize=8)

        if cutoff == 0:
            pylab.suptitle('%s' % (plottitle), size=11)
        else:
            pylab.suptitle('%s\nR < %d km' % (plottitle, cutoff), size=11)
        outfile = os.path.join(outdir, "gof-%s.png" % (gof_fileroot))
        print("==> Created GoF plot: %s" % (outfile))
        pylab.savefig(outfile, format="png",
                      transparent=False, dpi=plot_config.dpi)
        pylab.close()

def usage():
    """
    Prints usage information
    """
    print("usage: %s <title> <gof_fileroot> <indir> <outdir> <cutoff>" %
          (sys.argv[0]))

if __name__ == '__main__':
    if len(sys.argv) != 6:
        usage()
        sys.exit(1)

    MODE = 'rd50'
    PLOTTITLE = sys.argv[1]
    GOF_FILEROOT = sys.argv[2]
    INDIR = sys.argv[3]
    OUTDIR = sys.argv[4]
    CUTOFF = int(sys.argv[5])

    PLOTTER = PlotGoF()
    PLOTTER.plot(PLOTTITLE, GOF_FILEROOT, INDIR, OUTDIR, CUTOFF, MODE)
    sys.exit(0)
