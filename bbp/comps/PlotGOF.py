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

# Component extentions
COMP_EXT_RD50 = ['psa5n', 'psa5e', 'rotd50']
COMP_EXT_RD100 = ['rotd100', 'rotd50', 'ratio']

# Component titles
COMP_TITLE_RD50 = ['PSA North 5%', 'PSA East 5%', 'RotD50']
COMP_TITLE_RD100 = ['RotD100', 'RotD50', 'RotD100/RotD50']

# Component subplot locations
COMP_OFFSET = [312, 313, 311]

# Colorsets for plots
COLORSETS = {"single": 0,
             "combined": 1}
BIAS_COLORS = ["red", "black"]
BIAS_LH_COLORS = ["cyan", (0.8, 0.8, 1)]
CONF_LH_COLORS = ["yellow", (1, 0.3, 0.9)]

class PlotGoF(object):
    def __init__(self):
        return

    def read_data(self, datafile, min_period):
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
                (period <= MAX_PERIOD)):
                data[0].append(float(tmp[0]))
                data[1].append(float(tmp[1]))
        # Close file
        in_file.close()
        # Return data
        return data

    def multi_plot(self, plottitle, gof_fileroot, indir, outdir,
                   legends, num_stations):
        """
        Creates several GOF plots using the files specified in the
        gof_fileroot array
        """
        # Pick components and labels
        xtick_loc = XTICK_LOC_0_01
        xtick_label = XTICK_LABEL_0_01

        # Initialize data arrays
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
        for compnum in xrange(0, len(gof_fileroot)):
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
            for i in xrange(0, len(bias[compnum])):
                bias_l[compnum].append(bias[compnum][i] - sigma0[compnum][i])
                bias_h[compnum].append(bias[compnum][i] + sigma0[compnum][i])
                conf_l[compnum].append(m90[compnum][i])
                conf_h[compnum].append(p90[compnum][i])

        num_periods = len(periods[0])
        for comp in periods:
            if len(comp) != num_periods:
                print("Number of data points unequal across components")
                return

        # Construct baseline
        baseline = [0.0 for _ in periods[0]]

        # Find max, min values
        min_x = min([min(comp) for comp in periods])
        max_x = max([max(comp) for comp in periods])
        min_y = MIN_Y_AXIS
        max_y = MAX_Y_AXIS

        # Start plots
        num_plots = len(gof_fileroot)
        # Make 1 column of "num_plots" rows
        fig, axs = pylab.plt.subplots(num_plots, 1)
        # Set plot dims
        fig.set_size_inches(6, 10)

        #subplots_adjust(left=0.125)
        #subplots_adjust(right=0.9)
        #subplots_adjust(top=0.85)
        fig.subplots_adjust(hspace=0.4)
        fig.subplots_adjust(wspace=0.5)

        # Add subplots in a list
        subfigs = []
        for idx in range(0, num_plots):
            subfigs.append(axs[idx])

        # Now walk through each subfig
        for (subfig, subplot_title, cur_period,
             cur_bias, cur_bias_h, cur_bias_l,
             cur_conf_h, cur_conf_l) in zip(subfigs,
                                            legends,
                                            periods,
                                            bias,
                                            bias_h,
                                            bias_l,
                                            conf_h,
                                            conf_l):
            subfig.set_xlim(min_x, max_x)
            subfig.set_ylim(min_y, max_y)
            subfig.set_title("%s" % subplot_title, size=10)
            subfig.plot(cur_period, cur_bias, color='red', label='_nolegend_')
            subfig.fill_between(cur_period, cur_bias_h,
                                cur_bias_l, color='cyan',
                                label='_nolegend_')
            subfig.fill_between(cur_period, cur_conf_h,
                                cur_conf_l, color='yellow',
                                label='_nolegend_')
            subfig.plot(cur_period, baseline, color='grey',
                        label='_nolegend_')
            # Only put xlabel on bottom plot
            if legends.index(subplot_title) == len(gof_fileroot) - 1:
                subfig.set_xlabel("Period (sec)", size=8)
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
        outfile = os.path.join(outdir, "gof-%s.png" % (outfile))
        print("==> Created GoF plot: %s" % (outfile))
        fig.savefig(outfile, format="png",
                    transparent=False, dpi=plot_config.dpi)
        pylab.close()

    def plot(self, plottitle, gof_fileroot, indir, outdir,
             cutoff=0, min_period=0.01, mode=None, colorset=None):
        """
        Creates the GOF plot
        """

        # Pick components and labels
        start_comp = 0

        if min_period == 0.01:
            xtick_loc = XTICK_LOC_0_01
            xtick_label = XTICK_LABEL_0_01
        elif min_period == 0.1:
            xtick_loc = XTICK_LOC
            xtick_label = XTICK_LABEL
        else:
            raise bband_utils.ParameterError("invalid min_period: %f" %
                                             (min_period))

        if mode == "rd50" or mode == "rd50-single":
            comp_ext = COMP_EXT_RD50
            comp_title = COMP_TITLE_RD50
        elif mode == "rd100":
            comp_ext = COMP_EXT_RD100
            comp_title = COMP_TITLE_RD100
        else:
            raise bband_utils.ParameterError("plot mode %s unsupported" %
                                             (mode))
        #if mode == "rd50-single" or mode == "rd100":
        if mode == "rd50-single":
            # Only plot the last component
            # rd50 for rd50-single
            # rd100/rd50 ratio for rd100
            start_comp = 2

        # Select colorset
        if colorset is None:
            colorset_idx = 0
        else:
            colorset_idx = COLORSETS[colorset]

        # Set plot dims
        pylab.gcf().set_size_inches(6, 8)
        pylab.gcf().clf()

        pylab.subplots_adjust(left=0.125)
        pylab.subplots_adjust(right=0.9)
        pylab.subplots_adjust(top=0.85)
        pylab.subplots_adjust(hspace=0.4)
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

        for compnum in xrange(0, len(comp_ext)):
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
            for i in xrange(0, len(bias[compnum])):
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
        for i in xrange(0, len(period[0])):
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
        for compnum in xrange(start_comp, 3):
            comp = comp_ext[compnum]
            offset = COMP_OFFSET[compnum]

            pylab.subplot(offset)
            pylab.title(comp_title[compnum], size=10)
            pylab.plot(period[compnum], bias[compnum],
                       color=BIAS_COLORS[colorset_idx],
                       label='_nolegend_')
            pylab.fill_between(period[compnum], bias_h[compnum],
                               bias_l[compnum],
                               color=BIAS_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.fill_between(period[compnum], conf_h[compnum],
                               conf_l[compnum],
                               color=CONF_LH_COLORS[colorset_idx],
                               label='_nolegend_')
            pylab.plot(period[compnum], baseline, color='grey',
                       label='_nolegend_')
            pylab.xlim(min_x, max_x)
            if comp == 'ratio':
                pylab.ylim(MIN_Y_AXIS_RATIO, MAX_Y_AXIS_RATIO)
            else:
                pylab.ylim(min_y, max_y)
            if compnum == 1:
                pylab.xlabel("Period (sec)", size=8)
            pylab.ylabel("ln (data/model)", size=8)
            pylab.xscale('log')
            pylab.xticks(xtick_loc, xtick_label)
            pylab.tick_params(labelsize=8)

        if cutoff == 0:
            pylab.suptitle('%s' % (plottitle), size=12)
        else:
            pylab.suptitle('%s\nR < %d km' % (plottitle, cutoff), size=12)
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
    if (len(sys.argv) != 6):
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
