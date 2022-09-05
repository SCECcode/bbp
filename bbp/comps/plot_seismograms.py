#!/usr/bin/python
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

Broadband module to plot seismograms and overlay graphs
"""
from __future__ import division, print_function

# Import Python modules
import sys
import matplotlib as mpl
mpl.use("AGG")
import pylab

# Import plot config file
import plot_config

# S-wave velocity in km/s
S_VELOCITY = 4

def calculate_x_coords(ts1, rrup):
    """
    Calculated the min_x and max_x points using the timestamps and
    rrup as references
    """
    plot_mode = plot_config.plot_seismograms_mode
    max_ts = max(ts1)
    if plot_mode == 2:
        # Plot entire seismogram
        return (0, max_ts)
    if plot_mode == 1:
        # Plot first "duration" seconds
        return (0, plot_config.plot_seismograms_duration)
    if max_ts <= plot_config.plot_seismograms_duration:
        # Simulation is shorter than DURATION, plot everything
        return (0, max_ts)
    if rrup is None:
        # R_rup not available, plot first DURATION seconds
        return (0, min(max_ts, plot_config.plot_seismograms_duration))
    # R_rup provided, use it to calculate plot window
    min_x = float(rrup) / S_VELOCITY - 20.0
    if min_x < 0:
        min_x = 0
    max_x = min_x + plot_config.plot_seismograms_duration
    return (min_x, max_x)

def plot_seis(stat, filename, label, units, outfile, rrup=None):
    """
    Plots the seismogram for station stat, and outputs a png file outfile
    """
    ts1 = []
    ns1 = []
    ew1 = []
    ud1 = []

    cmt1 = ["", ""]

    # Read input file
    input_file = open(filename, 'r')
    for data in input_file:
        # Remove leading spaces
        data = data.strip()
        # Skip comments
        if data.startswith('#') or data.startswith('%'):
            if cmt1[0] == "":
                cmt1[0] = data
        else:
            tmp = []
            tmp = data.split()
            ts1.append(float(tmp[0]))
            ns1.append(float(tmp[1]))
            ew1.append(float(tmp[2]))
            ud1.append(float(tmp[3]))
    # Don't forget to close the file
    input_file.close()

    min_x, max_x = calculate_x_coords(ts1, rrup)
    min_horiz_y = 1.1 * min([min(ns1), min(ew1)])
    max_horiz_y = 1.1 * max([max(ns1), max(ew1)])
    min_vert_y = 1.1 * min(ud1)
    max_vert_y = 1.1 * max(ud1)

    # Set up ticks to match matplotlib 1.x style
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    pylab.clf()
    pylab.suptitle('Run %s, station %s' %
                   (label, stat), size=14)
    pylab.subplots_adjust(top=0.925)
    pylab.subplots_adjust(bottom=0.07)
    pylab.subplots_adjust(left=0.11)
    pylab.subplots_adjust(right=0.975)
    pylab.subplots_adjust(hspace=0.5)
    pylab.subplots_adjust(wspace=0.3)

    ax1 = pylab.subplot(311)
    pylab.plot(ts1, ns1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_horiz_y, max_horiz_y)
    if units == 'dis':
        pylab.ylabel("Displacement (cm)")
    elif units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")
    pylab.xlabel("Time (s)")
    ax1.set_title("N/S", fontsize="small")

    ax2 = pylab.subplot(312)
    pylab.plot(ts1, ew1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_horiz_y, max_horiz_y)
    if units == 'dis':
        pylab.ylabel("Displacement (cm)")
    elif units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")
    pylab.xlabel("Time (s)")
    ax2.set_title("E/W", fontsize="small")

    ax3 = pylab.subplot(313)
    pylab.plot(ts1, ud1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_vert_y, max_vert_y)
    if units == 'dis':
        pylab.ylabel("Displacement (cm)")
    elif units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")
    pylab.xlabel("Time (s)")
    ax3.set_title("U/D", fontsize="small")

    pylab.gcf().set_size_inches(6, 7)
    #pylab.tight_layout()
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

def read_seismogram_file(filename):
    """
    This function reads a seismogram and returns 4 lists with the
    horizontal components (ns and ew), vertical, and the timestamps
    """
    # Start empty
    ts = []
    ns = []
    ew = []
    ud = []

    # Read file
    seis_file = open(filename, 'r')
    for line in seis_file:
        # Remove leading spaces
        line = line.strip()
        # Skip comments
        if line.startswith('#') or line.startswith('%'):
            continue
        tmp = line.split()
        if len(tmp) < 4:
            print("Error reading seismogram in file %s" % (filename))
            sys.exit(1)
        ts.append(float(tmp[0]))
        ns.append(float(tmp[1]))
        ew.append(float(tmp[2]))
        ud.append(float(tmp[3]))
    # Close file
    seis_file.close()
    # All done
    return (ts, ns, ew, ud)

def plot_overlay(stat, obs_filename, comp_filename, obs_label, comp_label,
                 outfile, y_label="Velocity (cm/s)",
                 goflabel=None, gofdata=None):
    """
    This function plots observed and computed seismograms side by side
    for easy comparison
    """
    # Set up ticks to match matplotlib 1.x style
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    # Initialize variables
    textx = 0.53
    texty = 0.05
    fig = pylab.plt.figure()
    fig.clf()

    ts1, ns1, ew1, ud1 = read_seismogram_file(obs_filename)
    ts2, ns2, ew2, ud2 = read_seismogram_file(comp_filename)

    # Determine min and max X and Y for N/S/E/W/U/D for scaling
    min_x = 0
    max_x = min(max([max(ts1), max(ts2)]), 100)
    min_horiz_y = 1.1 * min([min(ns1), min(ns2), min(ew1), min(ew2)])
    max_horiz_y = 1.1 * max([max(ns1), max(ns2), max(ew1), max(ew2)])
    # Adjust so min and max are equal
    if abs(min_horiz_y) > abs(max_horiz_y):
        max_horiz_y = -1 * min_horiz_y
    else:
        min_horiz_y = -1 * max_horiz_y

    min_vert_y = 1.1 * min([min(ud1), min(ud2)])
    max_vert_y = 1.1 * max([max(ud1), max(ud2)])

    if abs(min_vert_y) > abs(max_vert_y):
        max_vert_y = -1 * min_vert_y
    else:
        min_vert_y = -1 * max_vert_y
    if goflabel is None or gofdata is None:
        fig.suptitle('%s vs %s, station %s' % (obs_label, comp_label, stat), size=14)
    else:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[0])
        fig.suptitle('%s vs %s, station %s (%s)' %
                     (obs_label, comp_label, stat, txt), size=14)
    fig.subplots_adjust(top=0.85)
    fig.subplots_adjust(left=0.075)
    fig.subplots_adjust(right=0.925)
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.3)

    ax = fig.add_subplot(231, title='%s, N/S' % obs_label)
    ax.plot(ts1, ns1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_ylabel(y_label)
    ax = fig.add_subplot(234, title='%s, N/S' % comp_label)
    ax.plot(ts2, ns2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_ylabel(y_label)

    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[2])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))

    #legend(prop=matplotlib.font_manager.FontProperties(size=10))

    ax = fig.add_subplot(232, title='%s, E/W' % obs_label)
    ax.plot(ts1, ew1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    #ylabel(y_label)
    ax = fig.add_subplot(235, title='%s, E/W' % comp_label)
    ax.plot(ts2, ew2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[1])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))
    #ylabel(y_label)
    #legend(prop=matplotlib.font_manager.FontProperties(size=10))

    ax = fig.add_subplot(233, title='%s, U/D' % obs_label)
    ax.plot(ts1, ud1, color='black', label=obs_label,
            lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_vert_y, max_vert_y)
    #ylabel(y_label)
    ax = fig.add_subplot(236, title='%s, U/D' % comp_label)
    ax.plot(ts2, ud2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_vert_y, max_vert_y)
    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[3])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))
    #ylabel(y_label)
    #legend(prop=matplotlib.font_manager.FontProperties(size=10))

    pylab.gcf().set_size_inches(10, 5)
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

def plot_overlay_with_arias(stat, obs_filename, comp_filename,
                            obs_arias_n_filename, obs_arias_e_filename,
                            obs_arias_z_filename, comp_arias_n_filename,
                            comp_arias_e_filename, comp_arias_z_filename,
                            obs_label, comp_label, outfile, rrup=None,
                            y_label="Velocity (cm/s)",
                            goflabel=None, gofdata=None):
    """
    This function plots observed and computed seismograms side by side
    for easy comparison
    """
    # Set up ticks to match matplotlib 1.x style
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    # Initialize variables
    textx = 0.53
    texty = 0.05
    fig = pylab.plt.figure()
    fig.clf()

    # Read all files
    (ts1, ns1, ew1, ud1) = read_seismogram_file(obs_filename)
    (ts2, ns2, ew2, ud2) = read_seismogram_file(comp_filename)
    ta1, _, _, an1 = read_seismogram_file(obs_arias_n_filename)
    ta1, _, _, ae1 = read_seismogram_file(obs_arias_e_filename)
    ta1, _, _, az1 = read_seismogram_file(obs_arias_z_filename)
    ta2, _, _, an2 = read_seismogram_file(comp_arias_n_filename)
    ta2, _, _, ae2 = read_seismogram_file(comp_arias_e_filename)
    ta2, _, _, az2 = read_seismogram_file(comp_arias_z_filename)

    # Determine min and max X and Y for N/S/E/W/U/D for scaling
    min_x = 0
    #max_x = min(max([max(ts1), max(ts2)]), 100)
    max_x = max([max(ts1), max(ts2)])
    min_horiz_y = 1.1 * min([min(ns1), min(ns2), min(ew1), min(ew2)])
    max_horiz_y = 1.1 * max([max(ns1), max(ns2), max(ew1), max(ew2)])

    # Adjust so min and max are equal
    if abs(min_horiz_y) > abs(max_horiz_y):
        max_horiz_y = -1 * min_horiz_y
    else:
        min_horiz_y = -1 * max_horiz_y

    min_vert_y = 1.1 * min([min(ud1), min(ud2)])
    max_vert_y = 1.1 * max([max(ud1), max(ud2)])

    if abs(min_vert_y) > abs(max_vert_y):
        max_vert_y = -1 * min_vert_y
    else:
        min_vert_y = -1 * max_vert_y
    # For arias plots, min=0, max=100%
    min_y_arias = 0
    max_y_arias = 100

    if goflabel is None or gofdata is None:
        fig.suptitle('%s vs %s, station %s' % (obs_label, comp_label, stat), size=14)
    else:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[0])
        fig.suptitle('%s vs %s, station %s (%s)' %
                     (obs_label, comp_label, stat, txt), size=14)
    fig.subplots_adjust(top=0.915)
    fig.subplots_adjust(left=0.075)
    fig.subplots_adjust(right=0.975)
    fig.subplots_adjust(bottom=0.07)
    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(wspace=0.2)

    ax = fig.add_subplot(321)
    ax.plot(ts1, ns1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_title("Observation N/S", fontsize='small')
    ax.set_ylabel(y_label)
    ax.set_xlabel("Time (s)")

    ax = fig.add_subplot(323)
    ax.plot(ts2, ns2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_title("Simulation N/S", fontsize='small')
    ax.set_ylabel(y_label)
    ax.set_xlabel("Time (s)")

    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[2])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))

    ax = fig.add_subplot(322)
    ax.plot(ts1, ew1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_title("Observation E/W", fontsize='small')
    ax.set_ylabel(y_label)
    ax.set_xlabel("Time (s)")

    ax = fig.add_subplot(324)
    ax.plot(ts2, ew2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_title("Simulation E/W", fontsize='small')
    ax.set_ylabel(y_label)
    ax.set_xlabel("Time (s)")

    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[1])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))

    ax = fig.add_subplot(325, title='N/S')
    ax.plot(ta1, an1, color='black', lw=plot_config.line_width)
    ax.plot(ta2, an2, color='red', lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y_arias, max_y_arias)
    ax.set_title("N/S", fontsize='small')
    ax.set_ylabel("Norm Arias Int (%)")
    ax.set_xlabel("Time (s)")

    ax = fig.add_subplot(326, title='E/W')
    ax.plot(ta1, ae1, color='black', lw=plot_config.line_width)
    ax.plot(ta2, ae2, color='red', lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y_arias, max_y_arias)
    ax.set_title("E/W", fontsize='small')
    ax.set_ylabel("Norm Arias Int (%)")
    ax.set_xlabel("Time (s)")

    pylab.gcf().set_size_inches(10, 7.5)
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

if __name__ == '__main__':
    plot_overlay(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                 sys.argv[5], sys.argv[6])
    #plot_seis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
