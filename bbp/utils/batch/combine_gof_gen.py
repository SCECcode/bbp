#!/usr/bin/env python
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

Python version of Ronnie Kamai Matlab scripts to generate a combined
bias plot. It collects information from the rd50 files for each
realization, groups the results by station (averaging) and then
generates a residuals file using the rd50 files from the recorded
data. This single residuals file uses the same resid2uncer_varN
program used in single bias plots to generate data for the combined
plot.
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

# Import Broadband modules
import bband_utils
from PlotGOF import PlotGoF
from gp_gof_cfg import GPGofCfg
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def summarize_rotd50(tmpdir, outdir, combined_resid_file,
                     comp_label, num_stations, num_realization,
                     codebase):
    """
    This function summarizes all rotd50 data and creates the combined
    rotd50 GOF plot
    """
    config = GPGofCfg()

    # Figure out where out binaries are
    if "BBP_DIR" in os.environ:
        install_root = os.path.normpath(os.environ["BBP_DIR"])
    else:
        raise bband_utils.ProcessingError("BBP_DIR is not set!")
    gp_bin_dir = os.path.join(install_root, "src", "gp", "bin")

    logfile = os.path.join(tmpdir, "log.txt")

    for comp in config.COMPS_PSA5:
        # Build paths and check lengths
        fileroot = os.path.join(tmpdir, "%s-%s-combined-rd50-%s" %
                                (codebase, comp_label, comp))
        bband_utils.check_path_lengths([combined_resid_file, fileroot],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s " % (os.path.join(gp_bin_dir, "resid2uncer_varN")) +
               "residfile=%s fileroot=%s " % (combined_resid_file, fileroot) +
               "comp=%s nstat=%d nper=63 " % (comp, num_stations) +
               " >> %s 2>&1" % (logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

    plottitle = ("Combined GOF Plot for %s\n%d Realizations - %s Method" %
                 (comp_label, num_realization, codebase.upper()))
    fileroot = "%s-%s-combined-rd50" % (codebase, comp_label)
    plotter = PlotGoF()
    plotter.plot(plottitle, fileroot, tmpdir, outdir,
                 cutoff=0, mode="rd50-single", colorset="combined")

    print("Stations used: %s" % (num_stations))

def combine_station_data(station, input_dir, temp_dir):
    """
    This function combines data for a given station across multiple
    realizations, writting a single output file in temp_dir
    """
    data = {}
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization)
        data_file = glob.glob("%s%s%s.%s.rd50" % (basedir, os.sep,
                                                  realization,
                                                  station))
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
            pieces = line.split()
            pieces = [float(piece) for piece in pieces]
            key = pieces[0]
            pieces = pieces[1:]
            if not key in data:
                # Key is new to dictionary
                empty_set = [[] for _ in pieces]
                data[key] = empty_set
            for idx, value in enumerate(pieces):
                data[key][idx].append(value)
        in_data.close()

    # Now, write the output file
    out_file = open((os.path.join(temp_dir, "%s.rd50" % (station))), 'w')
    keys = sorted(data.keys())
    for key in keys:
        out_file.write("%10.4f" % (key))
        for comp in data[key]:
            out_file.write(" %10.5e" % (numpy.mean(comp)))
        out_file.write("\n")

def combine_realizations_data(input_dir, temp_dir):
    """
    This function creates a single file averaging the rd50 files for
    each of the stations across all realizations
    """
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    one_realization = realizations[0]
    basedir = os.path.join(input_dir, one_realization)
    # Figure out what our stations are
    rd50_files = glob.glob("%s%s%s.*.rd50" % (basedir,
                                              os.sep,
                                              one_realization))
    rd50_files = [os.path.basename(each_file) for each_file in rd50_files]
    stations = [station.split(".")[1] for station in rd50_files]
    # Capture event_label
    bias_file = glob.glob("%s%s*.bias" % (basedir, os.sep))
    if len(bias_file) < 1:
        raise bband_utils.ProcessingError("Cannot find event label!")
    bias_file = bias_file[0]
    # Let's capture the event label
    event_label = os.path.basename(bias_file).split("-")[0]

    # Now walk through all realizations and combine stations data
    for station in stations:
        print("working on station: %s" % (station))
        combine_station_data(station, input_dir, temp_dir)

    return event_label, len(realizations), len(stations)

def create_resid_data_file(comp_label, input_indir, input_obsdir,
                           combined_file, temp_dir):
    """
    This function creates a file containing the combined residuals
    from the simulation data from all stations
    """
    # Copy header for first file, set logfile
    copy_header = 1
    logfile = os.path.join(temp_dir, "log.txt")

    # Figure out where out binaries are
    if "BBP_DIR" in os.environ:
        install_root = os.path.normpath(os.environ["BBP_DIR"])
    else:
        raise bband_utils.ProcessingError("BBP_DIR is not set!")
    gp_bin_dir = os.path.join(install_root, "src", "gp", "bin")

    # Get realizations
    realizations = sorted(os.listdir(input_indir))
    one_realization = realizations[0]
    basedir = os.path.join(input_indir, one_realization)

    # Get the station list
    a_statfile = glob.glob("%s%s*.stl" % (basedir, os.sep))
    if len(a_statfile) != 1:
        raise bband_utils.ProcessingError("Cannot get station list!")
    a_statfile = a_statfile[0]
    slo = StationList(a_statfile)
    site_list = slo.getStationList()

    # Get source file
    a_srcfile = glob.glob("%s%s*.src" % (basedir, os.sep))
    if len(a_srcfile) == 0:
        raise bband_utils.ProcessingError("Cannot get src file!")
    a_srcfile = a_srcfile[0]

    # Parse it!
    src_keys = bband_utils.parse_src_file(a_srcfile)

    # Get the obsdir
    realizations = sorted(os.listdir(input_obsdir))
    one_realization = realizations[0]
    basedir = os.path.join(input_obsdir, one_realization)
    obs_dir = glob.glob("%s%sobs_seis*" % (basedir, os.sep))
    if len(obs_dir) != 1:
        raise bband_utils.ProcessingError("Cannot get observation dir!")
    obs_dir = obs_dir[0]

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

        simfile1 = os.path.join(temp_dir, "%s.rd50" % (stat))
        datafile1 = os.path.join(obs_dir, "%s.rd50" % (stat))

        cmd = ("%s bbp_format=1 " %
               (os.path.join(gp_bin_dir, "gen_resid_tbl_3comp")) +
               "datafile1=%s simfile1=%s " % (datafile1, simfile1) +
               "comp1=psa5n comp2=psa5e comp3=rotd50 " +
               "eqname=%s mag=0.0 stat=%s lon=%.4f lat=%.4f " %
               (comp_label, stat, slon, slat) +
               "vs30=%d cd=%.2f " % (site.vs30, rrup) +
               "flo=%f fhi=%f " % (site.low_freq_corner,
                                   site.high_freq_corner) +
               "print_header=%d >> %s 2>> %s" %
               (copy_header, combined_file, logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

        if copy_header == 1:
            copy_header = 0


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
create_resid_data_file(COMP_LABEL, INPUT_INDIR, INPUT_OUTDIR,
                       COMBINED_FILE, TMPDIR)

summarize_rotd50(TMPDIR, OUTPUT_DIR,
                 COMBINED_FILE,
                 COMP_LABEL,
                 NUM_STAT,
                 NUM_REALIZATIONS,
                 OPTIONS.codebase)

print("All Done!")
# Clean-up, all done!
shutil.rmtree(TMPDIR)
