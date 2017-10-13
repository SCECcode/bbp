#!/usr/bin/python
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

Broadband module to plot seismograms and overlay graphs
"""
from __future__ import division, print_function

# Import Python modules
import sys
import matplotlib as mpl
mpl.use("AGG", warn=False)
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
    ver1 = []

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
            ver1.append(float(tmp[3]))
    # Don't forget to close the file
    input_file.close()

    min_x, max_x = calculate_x_coords(ts1, rrup)
    min_horiz_y = 1.1 * min([min(ns1), min(ew1)])
    max_horiz_y = 1.1 * max([max(ns1), max(ew1)])
    min_vert_y = 1.1 * min(ver1)
    max_vert_y = 1.1 * max(ver1)

    pylab.clf()
    pylab.suptitle('Seismograms for run %s, station %s' %
                   (label, stat), size=14)
    pylab.subplots_adjust(hspace=0.4)
    pylab.subplot(311, title='N/S')
    pylab.plot(ts1, ns1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_horiz_y, max_horiz_y)
    if units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")

    pylab.subplot(312, title='E/W')
    pylab.plot(ts1, ew1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_horiz_y, max_horiz_y)
    if units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")

    pylab.subplot(313, title='Ver')
    pylab.plot(ts1, ver1, lw=plot_config.line_width)
    pylab.xlim(min_x, max_x)
    pylab.ylim(min_vert_y, max_vert_y)
    if units == 'vel':
        pylab.ylabel("Velocity (cm/s)")
    elif units == 'acc':
        pylab.ylabel("Acceleration (cm/s/s)")

    pylab.gcf().set_size_inches(6, 7)
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
    ver = []

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
        ver.append(float(tmp[3]))
    # Close file
    seis_file.close()
    # All done
    return (ts, ns, ew, ver)

def plot_overlay(stat, obs_filename, comp_filename, obs_label, comp_label,
                 outfile, y_label="Velocity (cm/s)",
                 goflabel=None, gofdata=None):
    """
    This function plots observed and computed seismograms side by side
    for easy comparison
    """
    # Initialize variables
    textx = 0.53
    texty = 0.05
    fig = pylab.plt.figure()
    fig.clf()

    ts1, ns1, ew1, ver1 = read_seismogram_file(obs_filename)
    ts2, ns2, ew2, ver2 = read_seismogram_file(comp_filename)

    # Determine min and max X and Y for N/S/E/W, and Ver, for scaling
    min_x = 0
    max_x = min(max([max(ts1), max(ts2)]), 100)
    min_horiz_y = 1.1 * min([min(ns1), min(ns2), min(ew1), min(ew2)])
    max_horiz_y = 1.1 * max([max(ns1), max(ns2), max(ew1), max(ew2)])
    # Adjust so min and max are equal
    if abs(min_horiz_y) > abs(max_horiz_y):
        max_horiz_y = -1 * min_horiz_y
    else:
        min_horiz_y = -1 * max_horiz_y

    min_vert_y = 1.1 * min([min(ver1), min(ver2)])
    max_vert_y = 1.1 * max([max(ver1), max(ver2)])

    if abs(min_vert_y) > abs(max_vert_y):
        max_vert_y = -1 * min_vert_y
    else:
        min_vert_y = -1 * max_vert_y
    if goflabel is None or gofdata is None:
        fig.suptitle('%s vs. %s, station %s' % (obs_label, comp_label, stat), size=14)
    else:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[0])
        fig.suptitle('%s vs. %s, station %s (%s)' %
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
#       print "GOFLABEL, GOFDATA", goflabel, gofdata
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

    ax = fig.add_subplot(233, title='%s, ver' % obs_label)
    ax.plot(ts1, ver1, color='black', label=obs_label,
            lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_vert_y, max_vert_y)
    #ylabel(y_label)
    ax = fig.add_subplot(236, title='%s, ver' % comp_label)
    ax.plot(ts2, ver2, color='red', label=comp_label, lw=plot_config.line_width)
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
    # Initialize variables
    textx = 0.53
    texty = 0.05
    fig = pylab.plt.figure()
    fig.clf()

    # Read all files
    (ts1, ns1, ew1, ver1) = read_seismogram_file(obs_filename)
    (ts2, ns2, ew2, ver2) = read_seismogram_file(comp_filename)
    ta1, tmp1, tmp2, an1 = read_seismogram_file(obs_arias_n_filename)
    ta1, tmp1, tmp2, ae1 = read_seismogram_file(obs_arias_e_filename)
    ta1, tmp1, tmp2, az1 = read_seismogram_file(obs_arias_z_filename)
    ta2, tmp1, tmp2, an2 = read_seismogram_file(comp_arias_n_filename)
    ta2, tmp1, tmp2, ae2 = read_seismogram_file(comp_arias_e_filename)
    ta2, tmp1, tmp2, az2 = read_seismogram_file(comp_arias_z_filename)

    # Determine min and max X and Y for N/S/E/W, and Ver, for scaling
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

    min_vert_y = 1.1 * min([min(ver1), min(ver2)])
    max_vert_y = 1.1 * max([max(ver1), max(ver2)])

    if abs(min_vert_y) > abs(max_vert_y):
        max_vert_y = -1 * min_vert_y
    else:
        min_vert_y = -1 * max_vert_y
    # For arias plots, min=0, max=100%
    min_y_arias = 0
    max_y_arias = 100

    if goflabel is None or gofdata is None:
        fig.suptitle('%s vs. %s, station %s' % (obs_label, comp_label, stat), size=14)
    else:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[0])
        fig.suptitle('%s vs. %s, station %s (%s)' %
                     (obs_label, comp_label, stat, txt), size=14)
    fig.subplots_adjust(top=0.85)
    fig.subplots_adjust(left=0.075)
    fig.subplots_adjust(right=0.925)
    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(wspace=0.3)

    # FS: May 2013: for 3-comp plot below is #331
    ax = fig.add_subplot(321, title='%s, N/S' % obs_label)
    ax.plot(ts1, ns1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_ylabel(y_label)
    # FS: May 2013: for 3-comp plot below is #334
    ax = fig.add_subplot(323, title='%s, N/S' % comp_label)
    ax.plot(ts2, ns2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    ax.set_ylabel(y_label)
#       print "GOFLABEL, GOFDATA", goflabel, gofdata
    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[2])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))

    #legend(prop=matplotlib.font_manager.FontProperties(size=10))
    # FS: May 2013: for 3-comp plot below is #332
    ax = fig.add_subplot(322, title='%s, E/W' % obs_label)
    ax.plot(ts1, ew1, color='black', label=obs_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    #ylabel(y_label)

    # FS: May 2013: for 3-comp plot below is #335
    ax = fig.add_subplot(324, title='%s, E/W' % comp_label)
    ax.plot(ts2, ew2, color='red', label=comp_label, lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_horiz_y, max_horiz_y)
    if goflabel is not None and gofdata is not None:
        txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[1])
        ax.text(textx, texty, txt, transform=ax.transAxes,
                bbox=dict(facecolor='red', alpha=0.5))
    #ylabel(y_label)
    #legend(prop=matplotlib.font_manager.FontProperties(size=10))

    # FS: May 2013: Code commented out to remove vertical component
    # ax = fig.add_subplot(333, title='%s, ver' % obs_label)
    # ax.plot(ts1, ver1, color='black', label=obs_label,
    #         lw=plot_config.line_width)
    # ax.set_xlim(min_x, max_x)
    # ax.set_ylim(min_vert_y, max_vert_y)
    # #ylabel(y_label)
    # ax = fig.add_subplot(336, title='%s, ver' % comp_label)
    # ax.plot(ts2, ver2, color='red', label=comp_label,
    #         lw=plot_config.line_width)
    # ax.set_xlim(min_x, max_x)
    # ax.set_ylim(min_vert_y, max_vert_y)
    # if goflabel is not None and gofdata is not None:
    #     txt = '$%s_{%s}$=%.1f %%' % (goflabel[0], goflabel[1], gofdata[3])
    #     ax.text(textx, texty, txt, transform=ax.transAxes,
    #             bbox=dict(facecolor='red', alpha=0.5))

    #Ylabel(y_label)
    #legend(prop=matplotlib.font_manager.FontProperties(size=10))
    # FS: May 2013: for 3-comp plot below is #337
    ax = fig.add_subplot(325, title='Arias N/S')
    ax.plot(ta1, an1, color='black', lw=plot_config.line_width)
    ax.plot(ta2, an2, color='red', lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y_arias, max_y_arias)
    ax.set_ylabel("Norm Arias Int (%)")

    # FS: May 2013: for 3-comp plot below is #338
    ax = fig.add_subplot(326, title='Arias E/W')
    ax.plot(ta1, ae1, color='black', lw=plot_config.line_width)
    ax.plot(ta2, ae2, color='red', lw=plot_config.line_width)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y_arias, max_y_arias)

    # FS: May 2013: Code commented out to remove vertical component
    # ax = fig.add_subplot(339, title='Arias ver')
    # ax.plot(ta1, az1, color='black', lw=plot_config.line_width)
    # ax.plot(ta2, az2, color='red', lw=plot_config.line_width)
    # ax.set_xlim(min_x, max_x)
    # ax.set_ylim(min_y_arias, max_y_arias)

    pylab.gcf().set_size_inches(10, 7.5)
    pylab.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

if __name__ == '__main__':
    plot_overlay(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                 sys.argv[5], sys.argv[6])
    #plot_seis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
