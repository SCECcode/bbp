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

Program to merge a validation simulation with multiple segments with
several realizations
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import math
import random
import shutil
import argparse
import tempfile

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList
import validation_cfg

# Import BBP workflow modules
import plot_srf
#import plot_rotd50
from plot_seis import PlotSeis
from rotd50 import RotD50
from obs_seismograms import ObsSeismograms
from gen_plots import GenPlots
from gp_gof import GPGof

def post_process(station_list, src_files,
                 merged_outdir,
                 realization, val_obj):
    """
    Run the standard BBP post-processing tasks
    """
    # Plot seismograms
    plotter = PlotSeis(os.path.basename(station_list),
                       os.path.basename(src_files[0]),
                       True, True, int(realization))
    plotter.run()
    # RotD50
    process = RotD50(os.path.basename(station_list),
                     int(realization))
    process.run()
    process = ObsSeismograms(os.path.basename(station_list),
                             val_obj.get_obs_path(),
                             val_obj.get_obs_format(),
                             val_obj.get_obs_corrections(),
                             int(realization))
    process.run()
    process = GenPlots(os.path.basename(station_list),
                       val_obj.get_obs_path(),
                       'acc',
                       val_obj.get_validation_name(),
                       int(realization))
    process.run()
    process = GPGof(os.path.basename(src_files[0]),
                    os.path.basename(station_list),
                    val_obj.get_magnitude(),
                    val_obj.get_validation_name(),
                    val_obj.get_cutoff(),
                    False,
                    int(realization))
    process.run()

def plot_srf_file(merged_indir, merged_tmpdir, merged_outdir,
                  bbp_install, srf_files, val_event):
    """
    Creates the multi-segment SRF plot
    """
    # Save current directory
    old_cwd = os.getcwd()
    os.chdir(merged_tmpdir)

    for srf_file in srf_files:
        srf_file = os.path.basename(srf_file)
        # Write slip file
        srfbase = srf_file[0:srf_file.find(".srf")]
        slipfile = "%s.slip" % (srfbase)
        cmd = ("%s calc_xy=0 type=slip nseg=-1 < %s > %s" %
               (os.path.join(bbp_install.A_GP_BIN_DIR, "srf2xyz"),
                os.path.join(merged_indir, srf_file),
                slipfile))
        bband_utils.runprog(cmd)

        # Write tinit file
        tinitfile = "%s.tinit" % (srfbase)
        cmd = ("%s calc_xy=0 type=tinit nseg=-1 < %s > %s" %
               (os.path.join(bbp_install.A_GP_BIN_DIR, "srf2xyz"),
                os.path.join(merged_indir, srf_file),
                tinitfile))
        bband_utils.runprog(cmd)

    plottitle = 'Rupture Model for %s' % (val_event)
    plot_srf.plot_multi_segment(plottitle,
                                srf_files,
                                merged_outdir)

    # Restore directory
    os.chdir(old_cwd)

def add_bbp_seismograms(input_files, output_file):
    """
    Add all input files and write output_file with the combined data
    """
    # Start empty
    headers = []
    times = []
    ns_comp = []
    ew_comp = []
    ud_comp = []

    # Read first file
    input_file = input_files[0]
    input_files = input_files[1:]
    i_file = open(input_file, 'r')
    for line in i_file:
        line = line.strip()
        # Empty lines
        if not line:
            continue
        if line.startswith('#') or line.startswith('%'):
            headers.append(line)
            continue
        pieces = line.split()
        pieces = [float(piece) for piece in pieces]
        times.append(pieces[0])
        ns_comp.append(pieces[1])
        ew_comp.append(pieces[2])
        ud_comp.append(pieces[3])
    i_file.close()

    # Now add other files
    for input_file in input_files:
        index = 0
        i_file = open(input_file, 'r')
        for line in i_file:
            line = line.strip()
            # Empty lines
            if not line:
                continue
            if line.startswith('#') or line.startswith('%'):
                continue
            pieces = line.split()
            pieces = [float(piece) for piece in pieces]
            if index > len(times):
                print("[ERROR]: File size mismatch!")
                sys.exit(1)
            ns_comp[index] = ns_comp[index] + pieces[1]
            ew_comp[index] = ew_comp[index] + pieces[2]
            ud_comp[index] = ud_comp[index] + pieces[3]
            index = index + 1
        i_file.close()

    # Finally write output file
    o_file = open(output_file, 'w')
    for header_line in headers:
        o_file.write("%s\n" % (header_line))
    for time, ns_val, ew_val, ud_val in zip(times, ns_comp,
                                            ew_comp, ud_comp):
        o_file.write("%5.7f   %5.9e   %5.9e    %5.9e\n" %
                     (time, ns_val, ew_val, ud_val))
    o_file.close()

def merge_seismograms(input_sims, station_list,
                      merged_outdir, realization):
    """
    Adds seismograms from multiple simulations, creating a set
    of merged seismograms
    """
    # Load station list
    slo = StationList(station_list)
    site_list = slo.getStationList()

    # Merge each station
    for station in site_list:
        print("==> Merging station: %s" % (station.scode))
        # Merge both velocity and acceleration
        for file_type in ['vel', 'acc']:
            input_files = []
            for input_sim in input_sims:
                input_dir = os.path.join(input_sim, "Sims",
                                         "outdata", realization)
                input_file = os.path.join(input_dir,
                                          "%s.%s.%s.bbp" %
                                          (realization,
                                           station.scode,
                                           file_type))
                input_files.append(input_file)
                output_file = os.path.join(merged_outdir,
                                           "%s.%s.%s.bbp" %
                                           (realization,
                                            station.scode,
                                            file_type))
                add_bbp_seismograms(input_files, output_file)

def copy_indata_files(input_sims, merged_indir,
                      merged_tmpdir, realization):
    """
    Copies all needed indata files to the new indata directory
    """
    src_files = []
    srf_files = []
    for input_sim in input_sims:
        input_dir = os.path.join(input_sim, "Sims", "indata", realization)
        # SRC files
        src_file = glob.glob("%s/*.src" % (input_dir))
        if len(src_file) != 1:
            print("[ERROR]: Can't find single SRC file in %s" % (input_dir))
            sys.exit(1)
        src_file = src_file[0]
        src_files.append(src_file)
        # SRF files
        in_srf_file = "%s.srf" % (os.path.splitext(src_file)[0])
        tmp_srf_file = os.path.join(merged_tmpdir,
                                    os.path.basename(in_srf_file))
        srf_files.append(tmp_srf_file)
        # Copy files
        shutil.copy2(src_file, merged_indir)
        shutil.copy2(in_srf_file, merged_indir)
        shutil.copy2(in_srf_file, merged_tmpdir)

    # Now copy station list
    input_dir = os.path.join(input_sims[0], "Sims", "indata", realization)
    station_file = glob.glob("%s/*.stl" % (input_dir))
    if len(station_file) != 1:
        print("[ERROR]: Can't find station list file in %s" % (input_dir))
        sys.exit(1)
    station_file = station_file[0]
    station_list = station_file
    shutil.copy2(station_file, merged_indir)

    # Copy corrections file if needed
    correction_file = glob.glob("%s/*corrections.txt" % (input_dir))
    if len(correction_file) == 1:
        correction_file = correction_file[0]
        shutil.copy2(correction_file, merged_indir)

    return src_files, srf_files, station_list

def main():
    """
    Parse command line options and create the needed files/directories
    """
    prog_base = os.path.basename(sys.argv[0])
    usage = "usage: %s [options]" % (prog_base)
    parser = argparse.ArgumentParser(description="Merges a number of "
                                     "validation segments into a single "
                                     "multi-segment simulation")
    parser.add_argument("-d", "--dir", dest="merged_simdir",
                        help="Merged simulation directory")
    parser.add_argument("-e", "--event", dest="event",
                        help="validation event")
    parser.add_argument('input_segments', nargs='*',
                        help="top-level directories for all segments")

    args = parser.parse_args()

    # Input segments
    if len(args.input_segments) < 2:
        print("Please provide at least two simulations to merge!")
        sys.exit(-1)
    input_segments = args.input_segments

    # Event name
    if args.event is None:
        print("[ERROR]: Please provide an event name!")
        sys.exit(-1)
    event = args.event

    # Check for the simulation directory
    merged_simdir = args.merged_simdir
    if merged_simdir is None:
        print("Please provide a merged simulation directory!")
        sys.exit(1)
    merged_simdir = os.path.abspath(merged_simdir)
    if os.path.exists(merged_simdir):
        print("Merged simulation directory exists: %s" %
              (merged_simdir))
        opt = raw_input("Do you want to delete its contents (y/n)? ")
        if opt.lower() != "y":
            print("Please provide another simulation directory!")
            sys.exit(1)
        opt = raw_input("ARE YOU SURE (y/n)? ")
        if opt.lower() != "y":
            print("Please provide another simulation directory!")
            sys.exit(1)
        # Delete existing directory (we already asked the user twice!!!)
        shutil.rmtree(merged_simdir)

    # Create merged simulation directories
    os.makedirs(merged_simdir)
    merged_indir = os.path.join(merged_simdir, "Sims", "indata")
    merged_outdir = os.path.join(merged_simdir, "Sims", "outdata")
    merged_tmpdir = os.path.join(merged_simdir, "Sims", "tmpdata")
    merged_logdir = os.path.join(merged_simdir, "Sims", "logs")
    xmldir = os.path.join(merged_simdir, "Xml")
    srcdir = os.path.join(merged_simdir, "Src")
    for mdir in [merged_indir, merged_outdir, merged_tmpdir,
                 merged_logdir, xmldir, srcdir]:
        os.makedirs(mdir)

    # Setup BBP installation
    os.environ["BBP_DATA_DIR"] = os.path.join(merged_simdir, "Sims")
    bbp_install = InstallCfg.getInstance()

    # Make sure event exists
    event_names = validation_cfg.VE_EVENTS.get_all_names()
    event_names_lc = [e.lower() for e in event_names]
    if not event.lower() in event_names_lc:
        print("[ERROR]: Unknown event: %s" % (event))
        print("[ERROR]: Available events: ", event_names)
        shutil.rmtree(merged_simdir)
        sys.exit(1)
    val_event = event
    val_obj = validation_cfg.VE_EVENTS.get_event_by_print_name(event)

    # Make sure directories exist
    input_sims = []
    for simdir in input_segments:
        simdir = os.path.abspath(simdir)
        if not os.path.exists(simdir):
            print("Simulation directory %s does not exist!" % (simdir))
            sys.exit(1)
        input_sims.append(simdir)

    # Figure out the realizations
    input_sim = os.path.join(input_sims[0], "Sims", "outdata")
    realizations = glob.glob("%s%s*" % (input_sim, os.path.sep))
    num_realizations = len(realizations)
    realizations = [os.path.basename(item) for item in realizations]

    # Merge each realization
    print(" BBP Multi-segment merging tool ".center(80, "-"))
    print("=> Merging %d realizations..." % (num_realizations))
    for realization in realizations:
        print("==> Realization %s" % (realization))

        # Create directories
        merged_realization_indir = os.path.join(merged_indir,
                                                realization)
        merged_realization_tmpdir = os.path.join(merged_tmpdir,
                                                 realization)
        merged_realization_outdir = os.path.join(merged_outdir,
                                                 realization)
        merged_realization_logdir = os.path.join(merged_logdir,
                                                 realization)
        for mdir in [merged_realization_indir,
                     merged_realization_outdir,
                     merged_realization_tmpdir,
                     merged_realization_logdir]:
            os.makedirs(mdir)

        # Copy indata files
        (src_files,
         srf_files,
         station_list) = copy_indata_files(input_sims,
                                           merged_realization_indir,
                                           merged_realization_tmpdir,
                                           realization)

        # Merge seismograms
        merge_seismograms(input_sims, station_list,
                          merged_realization_outdir,
                          realization)
        # Plot combined SRF plot
        plot_srf_file(merged_realization_indir,
                      merged_realization_tmpdir,
                      merged_realization_outdir,
                      bbp_install,
                      srf_files,
                      val_event)
        # Post-processing steps
        post_process(station_list, src_files,
                     merged_realization_outdir,
                     realization,
                     val_obj)
        print("==> Realization %s Completed" % (realization))
        print("%s" % ("*" * 80))

if __name__ == "__main__":
    main()
