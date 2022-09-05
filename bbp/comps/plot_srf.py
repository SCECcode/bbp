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

Plots slip distribution for a SRF
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import numpy as np
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pylab

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg
from plot_utils import get_srf_num_segments, get_srf_params

# Import plot config file
import plot_config

# Slip range in cm
SLIP_X_FACTOR = 20.0
SLIP_Y_FACTOR = 5.0

def read_xy_file(input_file, numx, numy):
    """
    Read in fault file
    """
    my_file = open(input_file)
    slips = my_file.readlines()
    my_file.close()

    data = np.arange(numx * numy, dtype=float).reshape(numy, numx)

    # Data is x-fast
    for y in range(0, numy):
        for x in range(0, numx):
            tokens = slips[y * (numx) + x].split()
            data[y][x] = tokens[2]

    return data

def plot_multi_srf_files(plottitle, srffiles, outdir):
    """
    Produces the multi-segment SRF plot
    """
    num_segments = len(srffiles)

    srf_params = []
    srf_dims = []
    srf_extents = []
    srf_slips = []
    srf_tinits = []

    for srffile in srffiles:
        # Get SRF parameters
        params = get_srf_params(srffile)
        dim_len = params["dim_len"]
        dim_wid = params["dim_wid"]
        fault_len = params["fault_len"]
        fault_width = params["fault_width"]
        dims = [dim_len, dim_wid]
        extents = [-(fault_len / 2), (fault_len / 2),
                   fault_width, 0.0]

        # Read in SRF slips
        slipfile = "%s.slip" % (os.path.splitext(srffile)[0])
        slips = read_xy_file(slipfile, dims[0], dims[1])

        # Read in SRF tinits
        tinitfile = "%s.tinit" % (os.path.splitext(srffile)[0])
        tinits = read_xy_file(tinitfile, dims[0], dims[1])

        # Find avg/max slip
        sumslip = 0.0
        minslip = 100000.0
        maxslip = 0.0
        for y in range(0, dims[1]):
            for x in range(0, dims[0]):
                if slips[y][x] > maxslip:
                    maxslip = slips[y][x]
                if slips[y][x] < minslip:
                    minslip = slips[y][x]
                sumslip = sumslip + slips[y][x]

        params["minslip"] = minslip
        params["maxslip"] = maxslip
        params["sumslip"] = sumslip

        # Add to our lists
        srf_params.append(params)
        srf_dims.append(dims)
        srf_extents.append(extents)
        srf_slips.append(slips)
        srf_tinits.append(tinits)

    # Calculate min, max, average slip
    avgslip = 0.0
    minslip = 100000.0
    maxslip = 0.0
    totalpts = 0.0
    for params, dims in zip(srf_params, srf_dims):
        avgslip = avgslip + params["sumslip"]
        totalpts = totalpts + (dims[0] * dims[1])
        minslip = min(minslip, params["minslip"])
        maxslip = max(maxslip, params["maxslip"])
    avgslip = avgslip / totalpts

    # Create subfigures
    fig, subfigs = pylab.plt.subplots(1, num_segments, sharey=True)
    # Set plot dims
    fig.set_size_inches(11, 4.2)
    # Set title
    fig.suptitle('%s\nMin/Avg/Max Slip = %d/%d/%d' % (plottitle,
                                                      int(minslip),
                                                      int(avgslip),
                                                      int(maxslip)), size=12)

    # Set up propotions, first we calculate what we need
    num_spaces = num_segments - 1
    between_space = 0.02
    total_space = 0.9 # Figure goes from 0.05 to 0.95
    usable_space = total_space - between_space * num_spaces
    total_len = 0.0
    for params in srf_params:
        total_len = total_len + params["fault_len"]
    ratios = []
    for params in srf_params:
        ratios.append(params["fault_len"] / total_len)

    # Now we apply these to the axes
    current_position = 0.05
    for subfig, ratio in zip(subfigs, ratios):
        current_len = usable_space * ratio
        subfig.set_position([current_position, 0.2,
                             current_len, 0.60])
        current_position = current_position + current_len + between_space

    # Setup slip color scale
    cmap = cm.hot_r
    d = int(maxslip / SLIP_X_FACTOR + 0.0)
    while SLIP_X_FACTOR * d < 0.9 * maxslip:
        d = d + 1
    colormin = 0.0
    colormax = float(SLIP_X_FACTOR * d)
    colorint = float(SLIP_Y_FACTOR * d)
    norm = mcolors.Normalize(vmin=colormin, vmax=colormax)

    for (subfig, params, dims,
         slips, tinits,
         extents) in zip(subfigs, srf_params, srf_dims,
                         srf_slips, srf_tinits, srf_extents):

        subfig.set_adjustable('box')

        # Plot slips
        im = subfig.imshow(slips, cmap=cmap, norm=norm, extent=extents,
                           interpolation='nearest')

        # Freeze the axis extents
        subfig.set_autoscale_on(False)

        # Set font size
        for tick in subfig.get_xticklabels():
            tick.set_fontsize(8)
        for tick in subfig.get_yticklabels():
            tick.set_fontsize(8)

        subfig.set_title(u"Azimuth = %d\u00b0" % (params["azimuth"]), size=8)
        subfig.set_xlabel("Along Strike (km)", size=8)
        if subfig is subfigs[0]:
            subfig.set_ylabel("Down Dip (km)", size=8)

        # Setup tinit contours
        mintinit = 100000.0
        maxtinit = 0.0
        for y in range(0, dims[1]):
            for x in range(0, dims[0]):
                if tinits[y][x] > maxtinit:
                    maxtinit = tinits[y][x]
                if tinits[y][x] < mintinit:
                    mintinit = tinits[y][x]

        contour_intervals = ((maxtinit - mintinit) /
                             plot_config.PLOT_SRF_DEFAULT_CONTOUR_INTERVALS)
        # Plot tinit contours
        subfig.contour(tinits,
                       pylab.linspace(mintinit, maxtinit,
                                      int(round(contour_intervals))),
                       origin='upper', extent=extents, colors='k')

    # Setup slip color scale
    colorbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])
    cb = fig.colorbar(im, cax=colorbar_ax, orientation='horizontal',
                      ticks=pylab.linspace(colormin, colormax,
                                           int(colormax/colorint) + 1))
    cb.set_label('Slip (cm)', fontsize=8)
    for tick in cb.ax.get_xticklabels():
        tick.set_fontsize(8)

    # Save plot to file
    srffile = os.path.splitext(os.path.basename(srffiles[0]))[0]
    if srffile.find("_seg") > 0:
        srffile = srffile[0:srffile.find("_seg")]
    outfile = os.path.join(outdir, "%s.png" % (srffile))
    print("Saving plot to %s" % (outfile))
    pylab.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)

def plot_multi_plot(num_segments, srf_params, srf_dims,
                    srf_extents, srf_slips, srf_tinits,
                    plottitle, srffile, outdir):
    """
    Create actual plot for multi-segments
    """
    # Calculate min, max, average slip
    avgslip = 0.0
    minslip = 100000.0
    maxslip = 0.0
    totalpts = 0.0
    for params, dims in zip(srf_params, srf_dims):
        avgslip = avgslip + params["sumslip"]
        totalpts = totalpts + (dims[0] * dims[1])
        minslip = min(minslip, params["minslip"])
        maxslip = max(maxslip, params["maxslip"])
    avgslip = avgslip / totalpts

    # Create subfigures
    fig, subfigs = pylab.plt.subplots(1, num_segments, sharey=True)
    # Make sure it is an array
    if num_segments == 1:
        subfigs = [subfigs]
    # Set plot dims
    fig.set_size_inches(11, 4)
    # Set title
    fig.suptitle('%s\nMin/Avg/Max Slip = %d/%d/%d' % (plottitle,
                                                      int(minslip),
                                                      int(avgslip),
                                                      int(maxslip)), size=12)

    # Set up propotions, first we calculate what we need
    num_spaces = num_segments - 1
    between_space = 0.02
    total_space = 0.8 # Figure goes from 0.1 to 0.9
    usable_space = total_space - between_space * num_spaces
    total_len = 0.0
    for params in srf_params:
        total_len = total_len + params["fault_len"]
    ratios = []
    for params in srf_params:
        ratios.append(params["fault_len"] / total_len)

    # Now we apply these to the axes
    current_position = 0.1
    for subfig, ratio in zip(subfigs, ratios):
        current_len = usable_space * ratio
        subfig.set_position([current_position, 0.2,
                             current_len, 0.60])
        current_position = current_position + current_len + between_space

    # Setup slip color scale
    cmap = cm.hot_r
    d = int(maxslip / SLIP_X_FACTOR + 0.0)
    while SLIP_X_FACTOR * d < 0.9 * maxslip:
        d = d + 1
    colormin = 0.0
    colormax = float(SLIP_X_FACTOR * d)
    colorint = float(SLIP_Y_FACTOR * d)
    norm = mcolors.Normalize(vmin=colormin, vmax=colormax)

    for (subfig, params, dims,
         slips, tinits,
         extents) in zip(subfigs, srf_params, srf_dims,
                         srf_slips, srf_tinits, srf_extents):

        subfig.set_adjustable('box')

        # Plot slips
        im = subfig.imshow(slips, cmap=cmap, norm=norm, extent=extents,
                           interpolation='nearest')

        # Freeze the axis extents
        subfig.set_autoscale_on(False)

        # Set font size
        for tick in subfig.get_xticklabels():
            tick.set_fontsize(8)
        for tick in subfig.get_yticklabels():
            tick.set_fontsize(8)

        subfig.set_title(u"Azimuth = %d\u00b0" % (params["azimuth"]), size=8)
        subfig.set_xlabel("Along Strike (km)", size=8)
        if subfig is subfigs[0]:
            subfig.set_ylabel("Down Dip (km)", size=8)

        # Setup tinit contours
        mintinit = 100000.0
        maxtinit = 0.0
        for y in range(0, dims[1]):
            for x in range(0, dims[0]):
                if tinits[y][x] > maxtinit:
                    maxtinit = tinits[y][x]
                if tinits[y][x] < mintinit:
                    mintinit = tinits[y][x]

        contour_intervals = ((maxtinit - mintinit) /
                             plot_config.PLOT_SRF_DEFAULT_CONTOUR_INTERVALS)
        if contour_intervals < 10:
            contour_intervals = 10
        # Plot tinit contours
        subfig.contour(tinits,
                       pylab.linspace(mintinit, maxtinit,
                                      int(round(contour_intervals))),
                       origin='upper', extent=extents, colors='k')

    # Setup slip color scale
    colorbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])
    cb = fig.colorbar(im, cax=colorbar_ax, orientation='horizontal',
                      ticks=pylab.linspace(colormin, colormax,
                                           int((colormax/colorint) + 1)))
    cb.set_label('Slip (cm)', fontsize=8)
    for tick in cb.ax.get_xticklabels():
        tick.set_fontsize(8)

    # Save plot to file
    outfile = os.path.join(outdir,
                           "%s.png" %
                           (os.path.splitext(srffile)[0]))
    print("Saving plot to %s" % (outfile))
    pylab.savefig(outfile, format="png",
                  transparent=False, dpi=plot_config.dpi)

def plot(plottitle, srffile, outdir):
    """
    Produce the SRF plot
    """
    srf_params = []
    srf_dims = []
    srf_extents = []
    srf_slips = []
    srf_tinits = []

    # Get number of segments
    num_segments = get_srf_num_segments(srffile)

    for seg in range(num_segments):
        # Get SRF parameters
        params = get_srf_params(srffile, seg)
        dim_len = params["dim_len"]
        dim_wid = params["dim_wid"]
        fault_len = params["fault_len"]
        fault_width = params["fault_width"]

        # Plot dimensions
        dims = [dim_len, dim_wid]
        extents = [-(fault_len / 2), (fault_len / 2),
                   fault_width, 0.0]

        # Read in SRF slips
        slipfile = "%s_seg%d.slip" % (os.path.splitext(srffile)[0],
                                      seg)
        slips = read_xy_file(slipfile, dims[0], dims[1])

        # Read in SRF tinits
        tinitfile = "%s_seg%d.tinit" % (os.path.splitext(srffile)[0],
                                        seg)
        tinits = read_xy_file(tinitfile, dims[0], dims[1])

        # Find avg/max slip
        sumslip = 0.0
        minslip = 100000.0
        maxslip = 0.0
        for y in range(0, dims[1]):
            for x in range(0, dims[0]):
                if slips[y][x] > maxslip:
                    maxslip = slips[y][x]
                if slips[y][x] < minslip:
                    minslip = slips[y][x]
                sumslip = sumslip + slips[y][x]

        params["minslip"] = minslip
        params["maxslip"] = maxslip
        params["sumslip"] = sumslip

        # Add to our lists
        srf_params.append(params)
        srf_dims.append(dims)
        srf_extents.append(extents)
        srf_slips.append(slips)
        srf_tinits.append(tinits)

    if num_segments > 1:
        plot_multi_plot(num_segments, srf_params, srf_dims,
                        srf_extents, srf_slips, srf_tinits,
                        plottitle, srffile, outdir)
        return

    # Simple case for 1 segment only, keep it as before
    dims = srf_dims[0]
    tinits = srf_tinits[0]
    slips = srf_slips[0]
    extents = srf_extents[0]

    # Calculate min, max, average slip
    avgslip = 0.0
    minslip = 100000.0
    maxslip = 0.0
    totalpts = 0.0
    for params, dims in zip(srf_params, srf_dims):
        avgslip = avgslip + params["sumslip"]
        totalpts = totalpts + (dims[0] * dims[1])
        minslip = min(minslip, params["minslip"])
        maxslip = max(maxslip, params["maxslip"])
    avgslip = avgslip / totalpts

    # Set plot dims
    pylab.gcf().set_size_inches(6, 8)
    pylab.gcf().clf()

    # Set title and adjust title y-position
    t = pylab.title("%s\nAvg/Max Slip = %d/%d" % (plottitle,
                                                  int(avgslip),
                                                  int(maxslip)), size=12)
    t.set_y(1.05)

    # Setup slip color scale
    cmap = cm.hot_r
    d = int(maxslip / SLIP_X_FACTOR + 0.0)
    while SLIP_X_FACTOR * d < 0.9 * maxslip:
        d = d + 1
    colormin = 0.0
    colormax = float(SLIP_X_FACTOR * d)
    colorint = float(SLIP_Y_FACTOR * d)
    norm = mcolors.Normalize(vmin=colormin, vmax=colormax)

    # Plot slips
    pylab.imshow(slips, cmap=cmap,
                 norm=norm, extent=extents,
                 interpolation='nearest')

    # Freeze the axis extents
    pylab.gca().set_autoscale_on(False)
    pylab.xlabel("Along Strike (km)", size=8)
    pylab.ylabel("Down Dip (km)", size=8)

    # Set font size
    for tick in pylab.gca().get_xticklabels():
        tick.set_fontsize(8)
    for tick in pylab.gca().get_yticklabels():
        tick.set_fontsize(8)

    # Setup slip color scale
    cb = pylab.colorbar(orientation='horizontal', shrink=0.5,
                        ticks=pylab.linspace(colormin, colormax,
                                             int(colormax/colorint) + 1))
    cb.set_label('Slip (cm)', fontsize=8)
    for tick in cb.ax.get_xticklabels():
        tick.set_fontsize(8)

    # Setup tinit contours
    mintinit = 100000.0
    maxtinit = 0.0
    for y in range(0, dims[1]):
        for x in range(0, dims[0]):
            if tinits[y][x] > maxtinit:
                maxtinit = tinits[y][x]
            if tinits[y][x] < mintinit:
                mintinit = tinits[y][x]

    contour_intervals = ((maxtinit - mintinit) /
                         plot_config.PLOT_SRF_DEFAULT_CONTOUR_INTERVALS)
    if contour_intervals < 10:
        contour_intervals = 10
    # Plot tinit contours
    pylab.contour(tinits,
                  pylab.linspace(mintinit, maxtinit,
                                 int(round(contour_intervals))),
                  origin='upper', extent=extents, colors='k')

    outfile = os.path.join(outdir,
                           "%s.png" %
                           (os.path.splitext(srffile)[0]))
    print("Saving plot to %s" % (outfile))
    pylab.savefig(outfile, format="png",
                  transparent=False, dpi=plot_config.dpi)

def run(r_srffile, sim_id=0):
    """
    Creates a SRF plot from an SRF file
    """
    install = InstallCfg.getInstance()

    a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
    a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
    srf2xyz_bin = os.path.join(install.A_GP_BIN_DIR, "srf2xyz")

    # Save current directory
    old_cwd = os.getcwd()
    os.chdir(a_tmpdir)

    # Get number of segments
    num_segments = get_srf_num_segments(r_srffile)
    srfbase = r_srffile[0:r_srffile.find(".srf")]

    # Write slip and tinit files for each segment
    for seg in range(num_segments):
        slipfile = "%s_seg%d.slip" % (srfbase, seg)
        cmd = ("%s calc_xy=0 type=slip nseg=%d < %s > %s" %
               (srf2xyz_bin, seg, r_srffile, slipfile))
        bband_utils.runprog(cmd)

        tinitfile = "%s_seg%d.tinit" % (srfbase, seg)
        cmd = ("%s calc_xy=0 type=tinit nseg=%d < %s > %s" %
                (srf2xyz_bin, seg, r_srffile, tinitfile))
        bband_utils.runprog(cmd)

    plottitle = 'Rupture Model for %s' % (r_srffile)
    plot(plottitle, r_srffile, a_outdir)

    os.chdir(old_cwd)

def usage():
    """
    Prints program usage to the user
    """
    print("usage: %s srf_file <sim_id>" %
          (sys.argv[0]))
    return

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    SRF_FILE = sys.argv[1]
    SIMID = sys.argv[2]

    run(SRF_FILE, SIMID)
