#!/usr/bin/python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module contains functions to plot rotD50 output files
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
mpl.use('AGG', warn=False)
import pylab

# Import plot config file
import plot_config

def plot_rd50(stat, rotd50_file1, rotd50_file2, label1, label2,
              outfile, lfreq=None, hfreq=None, quiet=False):
    """
    This function generates the a plot with three subplots showing
    psa5 and rotd50 comparisons between two files
    """
    if not quiet:
        print("Plotting %s vs %s" % (rotd50_file1, rotd50_file2))

    periods1 = []
    pseudo_aa_ns1 = []
    pseudo_aa_ew1 = []
    rd50_aa1 = []

    periods2 = []
    pseudo_aa_ns2 = []
    pseudo_aa_ew2 = []
    rd50_aa2 = []

    rd50file1 = open(rotd50_file1, 'r')
    for line in rd50file1:
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        periods1.append(float(pieces[0]))
        pseudo_aa_ns1.append(float(pieces[1]))
        pseudo_aa_ew1.append(float(pieces[2]))
        rd50_aa1.append(float(pieces[3]))
    rd50file1.close()

    if periods1 == []:
        print("Input file %s is missing data! Aborting..." %
              (rotd50_file1))
        sys.exit(1)

    if rotd50_file2 != "-":
        rd50file2 = open(rotd50_file2, 'r')
        for line in rd50file2:
            if line.startswith("#") or line.startswith("%"):
                continue
            pieces = line.split()
            periods2.append(float(pieces[0]))
            pseudo_aa_ns2.append(float(pieces[1]))
            pseudo_aa_ew2.append(float(pieces[2]))
            rd50_aa2.append(float(pieces[3]))
        rd50file2.close()

        if periods2 == []:
            print("Input file %s is missing data! Aborting..." %
                  (rotd50_file2))
            sys.exit(1)
    else:
        # Only have one file to plot
        periods2 = periods1
        pseudo_aa_ns2 = pseudo_aa_ns1
        pseudo_aa_ew2 = pseudo_aa_ew1

    # Convert min/max frequencies to periods
    if lfreq is None:
        pmax = None
    else:
        pmax = 1.0 / float(lfreq)

    if hfreq is None:
        pmin = None
    else:
        pmin = 1.0 / float(hfreq)

    # Start plot
    pylab.clf()
    min_x = min([min((periods1)), min(periods2)])
    max_x = max([max((periods1)), max(periods2)])
    min_horiz_y = min([min(pseudo_aa_ns1), min(pseudo_aa_ns2),
                       min(pseudo_aa_ew1), min(pseudo_aa_ew2)]) / 1.1
    max_horiz_y = 1.1 * max([max(pseudo_aa_ns1), max(pseudo_aa_ns2),
                             max(pseudo_aa_ew1), max(pseudo_aa_ew2)])
    min_vert_y = min([min(pseudo_aa_ns1), min(pseudo_aa_ns2),
                      min(pseudo_aa_ew1), min(pseudo_aa_ew2)]) / 1.1
    max_vert_y = 1.1 * max([max(pseudo_aa_ns1), max(pseudo_aa_ns2),
                            max(pseudo_aa_ew1), max(pseudo_aa_ew2)])

    if rotd50_file2 != "-":
        pylab.suptitle('PSA for station %s, %s vs %s' % (stat, label1, label2),
                       size=14)
    else:
        pylab.suptitle('PSA for station %s' % (stat), size=14)

    pylab.subplots_adjust(top=0.85)
    pylab.subplots_adjust(bottom=0.15)
    pylab.subplots_adjust(left=0.075)
    pylab.subplots_adjust(right=0.975)
    pylab.subplots_adjust(hspace=0.3)
    pylab.subplots_adjust(wspace=0.3)

    # First plot
    ax1 = pylab.subplot(131)
    pylab.plot(periods1, pseudo_aa_ns1, label=str(label1))
    pylab.xlim(min_x, max_x)
    pylab.xscale('log')
    pylab.ylim(min_horiz_y, max_horiz_y)
    pylab.ylabel("PSA (g)")
    ax1.set_title('%s, N/S' % (label1), fontsize='small')
    if rotd50_file2 != "-":
        pylab.plot(periods2, pseudo_aa_ns2, label=str(label2))
    pylab.xlabel('Period (s)')
    if pmin is not None:
        pylab.vlines(pmin, min_horiz_y, max_horiz_y,
                     color='violet', linestyles='--')
    if pmax is not None:
        pylab.vlines(pmax, min_horiz_y, max_horiz_y, color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))

    # Second plot
    ax2 = pylab.subplot(132)
    pylab.plot(periods1, pseudo_aa_ew1, label=str(label1))
    pylab.xlim(min_x, max_x)
    pylab.xscale('log')
    pylab.ylim(min_horiz_y, max_horiz_y)
    pylab.ylabel("PSA (g)")
    ax2.set_title('%s, E/W' % (label1), fontsize='small')
    if rotd50_file2 != "-":
        pylab.plot(periods2, pseudo_aa_ew2, label=str(label2))
    pylab.xlabel('Period (s)')
    if pmin is not None:
        pylab.vlines(pmin, min_horiz_y, max_horiz_y,
                     color='violet', linestyles='--')
    if pmax is not None:
        pylab.vlines(pmax, min_horiz_y, max_horiz_y, color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))

    # Third plot
    ax3 = pylab.subplot(133)
    pylab.plot(periods1, rd50_aa1, label=str(label1))
    pylab.xlim(min_x, max_x)
    pylab.xscale('log')
    pylab.ylim(min_vert_y, max_vert_y)
    pylab.ylabel("PSA (g)")
    ax3.set_title('%s, RotD50' % (label1), fontsize='small')
    if rotd50_file2 != "-":
        pylab.plot(periods2, rd50_aa2, label=str(label2))
    pylab.xlabel('Period (s)')
    if pmin is not None:
        pylab.vlines(pmin, min_horiz_y, max_horiz_y,
                     color='violet', linestyles='--')
    if pmax is not None:
        pylab.vlines(pmax, min_horiz_y, max_horiz_y, color='r', linestyles='--')
    pylab.legend(prop=mpl.font_manager.FontProperties(size=8))
    pylab.gcf().set_size_inches(10, 4)
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

if __name__ == '__main__':
    PROG_BASE = os.path.split(sys.argv[0])[1]
    plot_rd50(sys.argv[2], sys.argv[3], sys.argv[4],
              sys.argv[5], sys.argv[6], sys.argv[7],
              sys.argv[8], sys.argv[9])
