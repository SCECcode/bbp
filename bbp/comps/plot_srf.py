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

# Import plot config file
import plot_config

# Slip range in cm
SLIP_X_FACTOR = 20.0
SLIP_Y_FACTOR = 5.0

def read_srf(input_file, numx, numy):
    """
    Read in fault file
    """
    my_file = open(input_file)
    slips = my_file.readlines()
    my_file.close()

    data = np.arange(numx * numy, dtype=float).reshape(numy, numx)

    # Data is x-fast
    for y in xrange(0, numy):
        for x in xrange(0, numx):
            tokens = slips[y * (numx) + x].split()
            data[y][x] = tokens[2]

    return data

def get_srf_params(srf_file):
    """
    Reads fault_len, width, dlen, and dwid from the srd file
    """
    srf_params1 = None
    srf_params2 = None
    srf = open(srf_file, 'r')
    for line in srf:
        if line.startswith("PLANE"):
            # Found the plane line, the next one should have what we need
            srf_params1 = srf.next()
            srf_params2 = srf.next()
            break
    srf.close()
    if srf_params1 is None or srf_params2 is None:
        print("Cannot determine parameters from SRF file %s" %
              (srf_file))
        sys.exit(1)
    srf_params1 = srf_params1.strip()
    srf_params1 = srf_params1.split()
    srf_params2 = srf_params2.strip()
    srf_params2 = srf_params2.split()
    # Make sure we have the correct number of pieces
    if len(srf_params1) != 6 or len(srf_params2) != 5:
        print("Cannot parse params from SRF file %s" %
              (srf_file))
        sys.exit(1)

    # Pick the parameters that we need
    params = {}
    params["dim_len"] = int(srf_params1[2])
    params["dim_wid"] = int(srf_params1[3])
    params["fault_len"] = float(srf_params1[4])
    params["fault_width"] = float(srf_params1[5])
    params["azimuth"] = int(float(srf_params2[0]))

    return params

def plot_multi_segment(plottitle, srffiles, outdir):
    """
    Produces the multi-segment SRF plot
    """
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
        slips = read_srf(slipfile, dims[0], dims[1])

        # Read in SRF tinits
        tinitfile = "%s.tinit" % (os.path.splitext(srffile)[0])
        tinits = read_srf(tinitfile, dims[0], dims[1])

        # Find avg/max slip
        sumslip = 0.0
        minslip = 100000.0
        maxslip = 0.0
        for y in xrange(0, dims[1]):
            for x in xrange(0, dims[0]):
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
    fig, subfigs = pylab.plt.subplots(1, len(srffiles), sharey=True)
    # Set plot dims
    fig.set_size_inches(11, 4)
    # Set title
    fig.suptitle('%s\nMin/Avg/Max Slip = %d/%d/%d' % (plottitle,
                                                      int(minslip),
                                                      int(avgslip),
                                                      int(maxslip)), size=12)

    # Set up propotions, first we calculate what we need
    num_spaces = len(srffiles) - 1
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

        subfig.set_adjustable('box-forced')

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
        for y in xrange(0, dims[1]):
            for x in xrange(0, dims[0]):
                if tinits[y][x] > maxtinit:
                    maxtinit = tinits[y][x]
                if tinits[y][x] < mintinit:
                    mintinit = tinits[y][x]

        contour_intervals = ((maxtinit - mintinit) /
                             plot_config.PLOT_SRF_DEFAULT_CONTOUR_INTERVALS)
        # Plot tinit contours
        subfig.contour(tinits,
                       pylab.linspace(mintinit, maxtinit,
                                      round(contour_intervals)),
                       origin='upper', extent=extents, colors='k')

    # Setup slip color scale
    colorbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])
    cb = fig.colorbar(im, cax=colorbar_ax, orientation='horizontal',
                      ticks=pylab.linspace(colormin, colormax,
                                           (colormax/colorint) + 1))
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

def plot(plottitle, srffile, outdir):
    """
    Produce the SRF plot
    """
    # Get SRF parameters
    params = get_srf_params(srffile)
    dim_len = params["dim_len"]
    dim_wid = params["dim_wid"]
    fault_len = params["fault_len"]
    fault_width = params["fault_width"]

    # Plot dimensions
    dims = [dim_len, dim_wid]
    extents = [-(fault_len / 2), (fault_len / 2),
               fault_width, 0.0]

    # Read in SRF slips
    slipfile = "%s.slip" % (os.path.splitext(srffile)[0])
    slips = read_srf(slipfile, dims[0], dims[1])

    # Read in SRF tinits
    tinitfile = "%s.tinit" % (os.path.splitext(srffile)[0])
    tinits = read_srf(tinitfile, dims[0], dims[1])

    # Find avg/max slip
    avgslip = 0.0
    minslip = 100000.0
    maxslip = 0.0
    for y in xrange(0, dims[1]):
        for x in xrange(0, dims[0]):
            if slips[y][x] > maxslip:
                maxslip = slips[y][x]
            if slips[y][x] < minslip:
                minslip = slips[y][x]
            avgslip = avgslip + slips[y][x]
    avgslip = avgslip / (dims[0] * dims[1])

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
    pylab.imshow(slips, cmap=cmap, norm=norm, extent=extents,
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
                                             (colormax/colorint) + 1))
    cb.set_label('Slip (cm)', fontsize=8)
    for tick in cb.ax.get_xticklabels():
        tick.set_fontsize(8)

    # Setup tinit contours
    mintinit = 100000.0
    maxtinit = 0.0
    for y in xrange(0, dims[1]):
        for x in xrange(0, dims[0]):
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
                                 round(contour_intervals)),
                  origin='upper', extent=extents, colors='k')

    outfile = os.path.join(outdir,
                           "%s.png" %
                           (os.path.splitext(srffile)[0]))
    print("Saving plot to %s" % (outfile))
    pylab.savefig(outfile, format="png",
                  transparent=False, dpi=plot_config.dpi)
    return

def run(r_srffile, sim_id=0):
    """
    Creates a SRF plot from an SRF file
    """
    install = InstallCfg.getInstance()

    a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
    a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))

    # Save current directory
    old_cwd = os.getcwd()
    os.chdir(a_tmpdir)

    # Write slip file
    srfbase = r_srffile[0:r_srffile.find(".srf")]
    slipfile = "%s.slip" % (srfbase)
    cmd = ("%s/srf2xyz calc_xy=0 type=slip nseg=-1 < %s > %s" %
           (install.A_GP_BIN_DIR, r_srffile, slipfile))
    bband_utils.runprog(cmd)

    # Write tinit file
    tinitfile = "%s.tinit" % (srfbase)
    cmd = ("%s/srf2xyz calc_xy=0 type=tinit nseg=-1 < %s > %s" %
           (install.A_GP_BIN_DIR, r_srffile, tinitfile))
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
    sys.exit(0)
