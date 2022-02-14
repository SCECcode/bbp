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

Script to combine FAS data from multiple realizations and create
a combined FAS GOF plot.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import numpy as np
import shutil
import optparse
import tempfile

# Import Broadband modules
import bband_utils
from PlotGOF import PlotGoF
from fas_gof_cfg import FASGofCfg, resid2uncer_py
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

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

def rewrite_fas_eas_file(fas_input_file, fas_output_file):
    """
    Reads the fas_input_file, and writes its
    content back without the eas column so that
    it can be used by the GoF tools
    """
    input_file = open(fas_input_file, 'r')
    output_file = open(fas_output_file, 'w')
    output_file.write("# Freq(Hz)\t FAS H1 (cm/s)\t FAS H2 (cm/s)\t "
                   "Smoothed EAS (cm/s)\n")
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
        output_file.write("%2.7E\t%2.7E\t%2.7E\t%2.7E\n" %
                          (pieces[0], pieces[1], pieces[2], pieces[4]))

    input_file.close()
    output_file.close()

def summarize_fas(tmpdir, outdir, fas_residfile,
                  comp_label, num_stations, num_realization,
                  codebase):
    """
    This function summarizes all fas data and creates the combined fas
    GOF plot

    """
    config = FASGofCfg()
    min_cdst = 0
    max_cdst = 1e+15

    logfile = os.path.join(tmpdir, "log.txt")

    for comp in config.COMPS_FAS:
        # Build paths and check lengths
        fileroot = os.path.join(tmpdir, "%s-%s-combined-fas-%s" %
                                (codebase, comp_label, comp))
        bband_utils.check_path_lengths([fas_residfile, fileroot],
                                       bband_utils.GP_MAX_FILENAME)
        resid2uncer_py(fas_residfile, fileroot, comp,
                                num_stations, min_cdst, max_cdst)

    plottitle = ("Combined GOF Plot for %s\n%d Realizations - %s Method" %
                 (comp_label, num_realization, codebase.upper()))
    fileroot = "%s-%s-combined-fas" % (codebase, comp_label)
    plotter = PlotGoF()
    plotter.plot_fas_gof(plottitle, fileroot, tmpdir, outdir,
                         cutoff=0, colorset="combined")

    print("Stations used: %s" % (num_stations))

def combine_station_data(station, outdir, tmpdir):
    """
    This function combines data for a given station across multiple
    realizations, writting a single output file in tmpdir
    """
    data = {}
    # Get realizations
    realizations = sorted(os.listdir(outdir))
    for realization in realizations:
        basedir = os.path.join(outdir, realization, "FAS")
        data_file = glob.glob("%s%s%s.%s.smc8.smooth.fs.col" % (basedir, os.sep,
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
    out_file = open((os.path.join(tmpdir, "%s.fas" % (station))), 'w')
    keys = sorted(data.keys())
    for key in keys:
        out_file.write("%10.4f" % (key))
        for comp in data[key]:
            out_file.write(" %10.5e" % (np.mean(comp)))
        out_file.write("\n")

def combine_realizations_data(outdir, tmpdir):
    """
    This function creates a single file averaging the FAS files for
    each of the stations across all realizations
    """
    # Get realizations
    realizations = sorted(os.listdir(outdir))
    one_realization = realizations[0]
    basedir = os.path.join(outdir, one_realization)
    # Figure out what our stations are
    fas_files = glob.glob("%s%sFAS%sobs*.col" % (basedir,
                                                 os.sep,
                                                 os.sep))
    fas_files = [os.path.basename(each_file) for each_file in fas_files]
    stations = [station.split(".")[1] for station in fas_files]
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
        combine_station_data(station, outdir, tmpdir)

    return event_label, len(realizations), len(stations)

def create_resid_data_file(comp_label, input_indir, input_outdir,
                           combined_file, tmpdir):
    """
    This function creates a file containing the combined residuals
    from the simulation data from all stations
    """
    # Print FAS header for first file, set logfile
    print_header_fas = 1
    logfile = os.path.join(tmpdir, "log.txt")

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

    # Get the fasdir
    realizations = sorted(os.listdir(input_outdir))
    one_realization = realizations[0]
    a_outdir = os.path.join(input_outdir, one_realization)
    fasdir = os.path.join(input_outdir, one_realization, "FAS")

    # Go through all stations
    for site in site_list:
        slon = float(site.lon)
        slat = float(site.lat)
        stat = site.scode

        # Pick up DT from simulated file
        acc_file = "%s.%s.acc.bbp" % (one_realization, stat)
        input_syn_acc_file = os.path.join(a_outdir, acc_file)
        syn_dt = read_bbp_dt(input_syn_acc_file)
        max_syn_freq = 1.0 / (syn_dt * 2)
        if max_syn_freq < site.high_freq_corner:
            print("station %s: freq: %f, syn_dt: %f" %
                  (stat, site.high_freq_corner, max_syn_freq))

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

        sim_file_in = os.path.join(tmpdir, "%s.fas" % (stat))
        obs_file_in = os.path.join(fasdir,
                                   "obs.%s.smc8.smooth.fs.col" %
                                   (stat))
        sim_file_tmp = os.path.join(tmpdir, "tmp.fas.sim.txt")
        obs_file_tmp = os.path.join(tmpdir, "tmp.fas.obs.txt")
        rewrite_fas_eas_file(sim_file_in, sim_file_tmp)
        rewrite_fas_eas_file(obs_file_in, obs_file_tmp)

        gen_resid_bin = os.path.join(gp_bin_dir, "gen_resid_tbl_3comp")
        
        cmd = ("%s bbp_format=1 " % (gen_resid_bin) +
               "datafile1=%s simfile1=%s " % (obs_file_tmp,
                                              sim_file_tmp) +
               "comp1=fash1 comp2=fash2 comp3=seas " +
               "eqname=%s mag=%s stat=%s lon=%.4f lat=%.4f " %
               (comp_label, src_keys['magnitude'], stat, slon, slat) +
               "vs30=%d cd=%.2f " % (site.vs30, rrup) +
               "flo=%f fhi=%f " % (1.0 / min(site.high_freq_corner,
                                             max_syn_freq),
                                   1.0 / site.low_freq_corner) +
               "print_header=%d >> %s 2>> %s" %
               (print_header_fas, combined_file, logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

        if print_header_fas == 1:
            print_header_fas = 0


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
                             "bbp-fas-resid-combined.txt")

# Combine realizations' data
(COMP_LABEL,
 NUM_REALIZATIONS,
 NUM_STAT) = combine_realizations_data(INPUT_OUTDIR,
                                       TMPDIR)

# Create data files with both gmpe and simulation data
create_resid_data_file(COMP_LABEL, INPUT_INDIR, INPUT_OUTDIR,
                       COMBINED_FILE, TMPDIR)

summarize_fas(TMPDIR, OUTPUT_DIR,
              COMBINED_FILE,
              COMP_LABEL,
              NUM_STAT,
              NUM_REALIZATIONS,
              OPTIONS.codebase)

print("All Done!")
# Clean-up, all done!
shutil.rmtree(TMPDIR)
