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

This program creates a combined plot for the RZZ2015 metrics,
using data from multiple realizations.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import numpy
import shutil
import optparse
import tempfile
import matplotlib as mpl
if (mpl.get_backend() != 'agg'):
    mpl.use('Agg') # Disable use of Tk/X11
import pylab

# Import Broadband modules
import bband_utils
import plot_config
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def plot_param(data, subfig, param, title, ylabel):
    """
    Plots a subplot using data from the data dict for the
    specified param
    """
    # Set label for observation
    obs_label = "%s_obs" % (param)
    gmpe_label = "%s_gmpe" % (param)
    subfig.set_title(title, size=14)
    for station in data:
        dist = data[station]["dist"]
        vals = data[station][param]
        vals_mean = numpy.array(numpy.mean(vals))
        vals_min = min(vals)
        vals_max = max(vals)
        obs = data[station][obs_label]
        gmpe = data[station][gmpe_label]
        subfig.plot(dist, obs, color='black', marker='o')
        subfig.set_ylabel(ylabel, size=10)
        # Plotting of the GMPE data is disabled right now
        # subfig.plot(dist, gmpe, color='blue', marker='s', fillstyle='none')
        subfig.errorbar(dist, vals_mean,
                        yerr=[[vals_mean - vals_min], [vals_max - vals_mean]],
                        color='red', marker='^',
                        linestyle='None')
    # subfig.set_ylabel("Acceleration, g", size=10)
    subfig.set_xlabel("Distance (km)", size=10)
    subfig.grid(True)
    subfig.minorticks_on()

def plot_rzz2015_data(tmpdir, outdir, data,
                      comp_label, num_stations, num_realization,
                      comp, codebase):
    """
    This function summarizes all rzz2015 data and creates a 5 panel
    plot with the R1-R5 parameters.
    """

    # Create path for out plot
    out_plot = os.path.join(outdir, "rzz2015-%s-%s-combined-%s.png" %
                            (comp_label, codebase, comp))

    # Create fig
    fig, axs = pylab.plt.subplots(2, 3)
    fig.set_size_inches(17, 8.5)
    fig.suptitle("RZZ2015 - %s - %s - %s" %
                 (comp_label, codebase, comp), size=16)
    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(left=0.05)
    fig.subplots_adjust(right=0.98)
    # Hide 6th subplot as we only have 5 to show
    axs[-1, -1].axis('off')

    # R1
    subfig = axs[0][0]
    plot_param(data, subfig, "r1", "Arias Intensity", "Ia (g^2.s)")

    # R2
    subfig = axs[0][1]
    plot_param(data, subfig, "r2", "Duration", "D5-95 (s)")

    # R3
    subfig = axs[0][2]
    plot_param(data, subfig, "r3", "Ia/D5-95", "Ia/D5-95 (g^2)")

    # R4
    subfig = axs[1][0]
    plot_param(data, subfig, "r4", "Wmid", "Wmid (Hz)")

    # R5
    subfig = axs[1][1]
    plot_param(data, subfig, "r5", "w'", "w' (Hz)")

    # All done! Save plot!
    pylab.savefig(out_plot, format="png", dpi=plot_config.dpi)

def combine_station_data(station, input_dir, temp_dir):
    """
    This function combines data for a given station across multiple
    realizations, writting a single output file in temp_dir
    """
    data = []
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization,
                               "validations", "rzz2015")
        data_file = glob.glob("%s%s%s.rzz2015.*.txt" % (basedir,
                                                        os.sep,
                                                        realization))
        if len(data_file) != 1:
            raise bband_utils.ProcessingError("Data for station %s " %
                                              (station) +
                                              "not found for "
                                              "realization %s!" %
                                              (realization))
        data_file = data_file[0]
        in_data = open(data_file, 'r')
        for line in in_data:
            line = line.strip()
            # Skip comments
            if line.startswith("#"):
                continue
            pieces = line.split(",")
            cur_station = pieces[0]
            # Check if this is the station we want
            if cur_station == station:
                data.append(line)
        in_data.close()

    # Now, write the output file
    out_file = open((os.path.join(temp_dir, "%s.rzz2015" % (station))), 'w')
    for item in data:
        out_file.write("%s\n" % (item))
    out_file.close()

def combine_realizations_data(input_dir, temp_dir):
    """
    This function creates a single file averaging the rd100 files for
    each of the stations across all realizations
    """
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    one_realization = realizations[0]
    basedir = os.path.join(input_dir, one_realization,
                           "validations", "rzz2015")
    basedir_gmpe = os.path.join(input_dir, one_realization,
                                "validations", "rzz2015_gmpe")
    # Figure out what our stations are
    rzzgmpe_files = glob.glob("%s%s%s.rzz2015gmpe.*.png" % (basedir_gmpe,
                                                            os.sep,
                                                            one_realization))
    rzzgmpe_files = [os.path.basename(each_file) for each_file in rzzgmpe_files]
    stations = [station.split(".")[2] for station in rzzgmpe_files]
    # Capture event_label
    rzz_file = glob.glob("%s%s%s.rzz2015.*.txt" % (basedir, os.sep,
                                                   one_realization))
    if len(rzz_file) < 1:
        raise bband_utils.ProcessingError("Cannot find event label!")
    rzz_file = rzz_file[0]
    # Let's capture the event label
    event_label = os.path.basename(rzz_file).split(".")[2]

    # Now walk through all realizations and combine stations data
    for station in stations:
        print("working on station: %s" % (station))
        combine_station_data(station, input_dir, temp_dir)

    return event_label, len(realizations), len(stations)

def load_all_data(comp_label, input_indir, input_obsdir,
                  combined_file, temp_dir, component):
    """
    This function loads all data from each station file
    and creates the structures needed for plotting.
    """
    data = {}

    # Get realizations
    realizations = sorted(os.listdir(input_indir))
    one_realization = realizations[0]
    basedir = os.path.join(input_indir, one_realization)

    # Get the GMPE data for the RZZ2015 metrics
    base_outdir = os.path.join(input_obsdir, one_realization,
                               "validations", "rzz2015_gmpe")
    a_rzz2015_gmpe = glob.glob("%s%s%s.rzz2015gmpe.txt" % (base_outdir,
                                                           os.sep,
                                                           one_realization))
    a_rzz2015_gmpe = a_rzz2015_gmpe[0]
    # Get the station list
    a_statfile = glob.glob("%s%s*.stl" % (basedir, os.sep))
    if len(a_statfile) != 1:
        raise bband_utils.ProcessingError("Cannot get station list!")
    a_statfile = a_statfile[0]
    slo = StationList(a_statfile)
    site_list = slo.getStationList()

    # Get source file
    a_srcfile = glob.glob("%s%s*.src" % (basedir, os.sep))
    if len(a_srcfile) != 1:
        raise bband_utils.ProcessingError("Cannot get src file!")
    a_srcfile = a_srcfile[0]

    # Parse it!
    src_keys = bband_utils.parse_src_file(a_srcfile)

    # Go through all stations
    for site in site_list:
        slon = float(site.lon)
        slat = float(site.lat)
        stat = site.scode

        # Calculate Rrup
        origin = (src_keys['lon_top_center'],
                  src_keys['lat_top_center'])
        dims = (src_keys['fault_length'], src_keys['dlen'],
                src_keys['fault_width'], src_keys['dwid'],
                src_keys['depth_to_top'])
        mech = (src_keys['strike'], src_keys['dip'],
                src_keys['rake'])

        site_geom = [float(site.lon), float(site.lat), 0.0]
        (fault_trace1, up_seis_depth,
         low_seis_depth, ave_dip,
         dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
        _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                         fault_trace1,
                                                         up_seis_depth,
                                                         low_seis_depth,
                                                         ave_dip)

        # Read data for this station
        data_file = os.path.join(temp_dir, "%s.rzz2015" % (stat))

        data[stat] = {}
        data[stat]["dist"] = rrup
        data[stat]["r1"] = []
        data[stat]["r2"] = []
        data[stat]["r3"] = []
        data[stat]["r4"] = []
        data[stat]["r5"] = []
        data[stat]["r1_obs"] = None
        data[stat]["r2_obs"] = None
        data[stat]["r3_obs"] = None
        data[stat]["r4_obs"] = None
        data[stat]["r5_obs"] = None
        data[stat]["r1_gmpe"] = None
        data[stat]["r2_gmpe"] = None
        data[stat]["r3_gmpe"] = None
        data[stat]["r4_gmpe"] = None
        data[stat]["r5_gmpe"] = None

        in_file = open(data_file, 'r')
        for line in in_file:
            line = line.strip()
            if line.startswith("#"):
                # Skip comments
                continue
            pieces = line.split(",")
            comp = pieces[1].strip()
            # Check if we want this component
            if component != "both":
                if comp != component:
                    # Skip
                    continue
            # We want this data point
            pieces = pieces[2:]
            pieces = [float(piece) for piece in pieces]
            # Get observation values
            if data[stat]["r1_obs"] is None:
                data[stat]["r1_obs"] = pieces[6]
            if data[stat]["r2_obs"] is None:
                data[stat]["r2_obs"] = pieces[8]
            if data[stat]["r3_obs"] is None:
                data[stat]["r3_obs"] = pieces[10]
            if data[stat]["r4_obs"] is None:
                data[stat]["r4_obs"] = pieces[12]
            if data[stat]["r5_obs"] is None:
                data[stat]["r5_obs"] = pieces[14]
            # Get simulated data values
            data[stat]["r1"].append(pieces[7])
            data[stat]["r2"].append(pieces[9])
            data[stat]["r3"].append(pieces[11])
            data[stat]["r4"].append(pieces[13])
            data[stat]["r5"].append(pieces[15])
        in_file.close()

    gmpe_file = open(a_rzz2015_gmpe, 'r')
    for line in gmpe_file:
        line = line.strip()
        # Skip comments
        if line.startswith("#"):
            continue
        pieces = line.split(",")
        stat = pieces[0].strip()
        pieces = pieces[1:]
        pieces = [float(piece.strip()) for piece in pieces]
        data[stat]["r1_gmpe"] = pieces[2]
        data[stat]["r2_gmpe"] = pieces[3]
        data[stat]["r3_gmpe"] = pieces[2]/pieces[3]
        data[stat]["r4_gmpe"] = pieces[5]
        data[stat]["r5_gmpe"] = pieces[6]
    gmpe_file.close()

    # Return all data
    return data

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

PARSER = optparse.OptionParser()
PARSER.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
PARSER.add_option("-o", "--output_dir", dest="output_dir",
                  help="Output file")
PARSER.add_option("-c", "--codebase", dest="codebase",
                  help="Method used for the simulation")
PARSER.add_option("-p", "--component", dest="component",
                  help="Component (default: both 1 and 2)")
(OPTIONS, ARGS) = PARSER.parse_args()

if OPTIONS.input_dir is None:
    PARSER.error("Please specify the input directory!")
TOP_INPUT_DIR = OPTIONS.input_dir
if not os.path.isdir(TOP_INPUT_DIR):
    PARSER.error("Invalid input directory!")
if not "Sims" in os.listdir(TOP_INPUT_DIR):
    PARSER.error("Please provide the top-level simulation directory!\n"
                 "This is the directory given to the cluster script")
INPUT_OUTDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "outdata")
INPUT_TMPDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "tmpdata")
INPUT_INDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "indata")

if OPTIONS.output_dir is None:
    PARSER.error("error specify output directory!")
else:
    OUTPUT_DIR = OPTIONS.output_dir
    if not os.path.isdir(OUTPUT_DIR):
        PARSER.error("Invalid output directory!")

if OPTIONS.codebase is None:
    PARSER.error("Please specify codebase!")

if OPTIONS.component is None:
    COMP = "both"
else:
    if int(float(OPTIONS.component)) == 1:
        COMP = "001"
    elif int(float(OPTIONS.component)) == 2:
        COMP = "002"
    else:
        PARSER.error("Component must be 1 or 2")

# Create temp dir
TMPDIR = tempfile.mkdtemp(prefix="bbp-")
COMBINED_FILE = os.path.join(TMPDIR,
                             "bbp-rd50-resid-combined.txt")

# Combine realizations' data
(COMP_LABEL,
 NUM_REALIZATIONS,
 NUM_STAT) = combine_realizations_data(INPUT_OUTDIR,
                                       TMPDIR)

# Create data files with both gmpe and simulation data
DATA = load_all_data(COMP_LABEL,
                     INPUT_INDIR,
                     INPUT_OUTDIR,
                     COMBINED_FILE,
                     TMPDIR, COMP)

plot_rzz2015_data(TMPDIR, OUTPUT_DIR,
                  DATA,
                  COMP_LABEL,
                  NUM_STAT,
                  NUM_REALIZATIONS,
                  COMP,
                  OPTIONS.codebase)

print("All Done!")
# Clean-up, all done!
shutil.rmtree(TMPDIR)
