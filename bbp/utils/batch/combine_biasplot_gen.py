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

Python version of Ronnie Kamai Matlab scripts to generate a combined
bias plot. It compiles information from the rd50-resid.txt files from
all realizations, writes a single temporary file with all the data and
uses the same resid2uncer_varN program used in single bias plots to
generate data for the combined plot.
"""

# Import Python modules
import os
import glob
import shutil
import optparse
import tempfile

# Import Broadband modules
from PlotGOF import PlotGoF
from gp_gof_cfg import GPGofCfg
import bband_utils

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

        cmd = ("%s/resid2uncer_varN " % (gp_bin_dir) +
               "residfile=%s fileroot=%s " % (combined_resid_file, fileroot) +
               "comp=%s nstat=%d nper=63 " % (comp, num_stations) +
               " >> %s 2>&1" % (logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

    plottitle = ("Combined GOF Plot for %s\n%d Realizations\n%s Method" %
                 (comp_label, num_realization, codebase.upper()))
    fileroot = "%s-%s-combined-rd50" % (codebase, comp_label)
    plotter = PlotGoF()
    plotter.plot(plottitle, fileroot, tmpdir, outdir,
                 cutoff=0, mode="rd50-single", colorset="combined")

    print "Stations used: %s" % (num_stations)

def create_resid_data_file(input_dir, combined_file):
    """
    This function creates a file containing the combined residuals
    from the simulation data from all realizations
    """
    copy_header = True
    event_label = None
    num_stations = 0

    # Open output file, write header
    outfile = open(combined_file, 'w')

    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization)
        resid_file = glob.glob("%s%s*.rd50-resid.txt" % (basedir, os.sep))
        if len(resid_file) != 1:
            raise bband_utils.ProcessingError("Residuals file not found for "
                                              "realization %s!" % (realization))
        resid_file = resid_file[0]
        # Let's capture the event label
        if event_label is None:
            event_label = os.path.basename(resid_file).split("-")[0]
        input_file = open(resid_file, 'r')
        for line in input_file:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith("#") or line.startswith("%") or not line:
                continue
            if line.startswith("EQ"):
                # This is the header line, skip if already done
                if not copy_header:
                    continue
                # If not done, set flag
                copy_header = False
            if line.find("rotd50") > 0:
                num_stations = num_stations + 1
            outfile.write("%s\n" % (line))
        input_file.close()

    # All done!
    outfile.close()

    # Return event label
    return event_label, len(realizations), num_stations

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

# Create data files with both gmpe and simulation data
COMP_LABEL, NUM_REALIZATIONS, NUM_STAT = create_resid_data_file(INPUT_OUTDIR,
                                                                COMBINED_FILE)
summarize_rotd50(TMPDIR, OUTPUT_DIR,
                 COMBINED_FILE,
                 COMP_LABEL,
                 NUM_STAT,
                 NUM_REALIZATIONS,
                 OPTIONS.codebase)

print "All Done!"
# Clean-up, all done!
shutil.rmtree(TMPDIR)
