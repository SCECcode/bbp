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

Python version of Ronnie Kamai Matlab scripts to generate box plots
of the GMPE comparisons. Still uses one of Ronnies Matlab scripts to
actually generate the plots, but Python code is used to create
intermediate data files containing data from the GMPE runs.
"""

# Import Python modules
import os
import sys
import glob
import shutil
import optparse
import tempfile

# Import Broadband modules
import bband_utils

# --------------------------------------------------------------------------
# Some constants
# --------------------------------------------------------------------------
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "irikura", "song"]
PERIODS = [0.010, 0.011, 0.012, 0.013, 0.015, 0.017, 0.020, 0.022, 0.025,
           0.029, 0.032, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065,
           0.075, 0.085, 0.100, 0.110, 0.120, 0.130, 0.150, 0.170, 0.200,
           0.220, 0.240, 0.260, 0.280, 0.300, 0.350, 0.400, 0.450, 0.500,
           0.550, 0.600, 0.650, 0.750, 0.850, 1.000, 1.100, 1.200, 1.300,
           1.500, 1.700, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.500,
           4.000, 4.400, 5.000, 5.500, 6.000, 6.500, 7.500, 8.500, 10.000]
BOXPLOT = "gmpe_boxplot"
# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def create_gmpe_data_file(input_dir, gmpe_file):
    """
    This function creates a file containing the GMPE data
    """
    # Pick any realization (the GMPE data is the same for all of them)
    try:
        rel_dir = os.listdir(input_dir)[0]
        input_dir = os.path.join(input_dir, rel_dir)
        # Now get the obs_gmpe directory
        input_dir = glob.glob("%s%sgmpe_data*" % (input_dir, os.sep))[0]
    except IndexError:
        print "Unable to find a GMPE data directory!"
        print "Perhaps this simulation is still running..."
        sys.exit(1)

    # Open output file, write header
    outfile = open(gmpe_file, 'w')
    #outfile.write("GMPE, Magnitude, Mechanism, Distance, Station Name,")
    # Add header for the GMPE column
    outfile.write("0")
    for period in PERIODS:
        outfile.write(",%10.5f" % period)
    outfile.write("\n")

    # Get list of stations to process
    stations = sorted(glob.glob("%s%s*-gmpe.ri50" % (input_dir, os.sep)))
    for station in stations:
        st_id = os.path.basename(station)
        # Find just the numeric part of the station id
        st_id = st_id[st_id.find('-')+1:st_id.rfind('-')]
        # Start empty
        AS = []
        BA = []
        CB = []
        CY = []
        input_file = open(station, 'r')
        for line in input_file:
            line = line.strip()
            # Skip comments
            if line.startswith("#"):
                continue
            pieces = [float(item) for item in line.split()]
            AS.append(pieces[1])
            BA.append(pieces[2])
            CB.append(pieces[3])
            CY.append(pieces[4])
        # Read all values
        input_file.close()
        #outfile.write("AS08, %2.1f, %s, %d, %s, " %
        #              (MAGNITUDE, MECH, DIST, st_id))
        outfile.write("1")
        for item in AS:
            outfile.write(",%10.6f" % (item))
        outfile.write("\n")
        #outfile.write("BA08, %2.1f, %s, %d, %s, " %
        #              (MAGNITUDE, MECH, DIST, st_id))
        outfile.write("2")
        for item in BA:
            outfile.write(",%10.6f" % (item))
        outfile.write("\n")
        #outfile.write("CB08, %2.1f, %s, %d, %s, " %
        #              (MAGNITUDE, MECH, DIST, st_id))
        outfile.write("3")
        for item in CB:
            outfile.write(",%10.6f" % (item))
        outfile.write("\n")
        #outfile.write("CY08, %2.1f, %s, %d, %s, " %
        #              (MAGNITUDE, MECH, DIST, st_id))
        outfile.write("4")
        for item in CY:
            outfile.write(",%10.6f" % (item))
        outfile.write("\n")
    # All done, close output file
    outfile.close()

def create_sim_data_file(input_dir, sim_file):
    """
    This function creates a file containing the simulation data
    """
    # Open output file, write header
    outfile = open(sim_file, 'w')
    #outfile.write("realization, Magnitude, Mechanism,"
    #              "Distance, Station Name,")
    outfile.write("%10.5f" % (PERIODS[0]))
    for period in PERIODS[1:]:
        outfile.write(",%10.5f" % period)
    outfile.write("\n")

    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization)
        stations = sorted(glob.glob("%s%s*.rd50" % (basedir, os.sep)))
        for station in stations:
            st_id = os.path.basename(station)
            # Find just the numeric part of the station id
            st_id = st_id[st_id.find('-')+1:st_id.find(".rd50")]
            input_file = open(station, 'r')
            # Write first part of line
            #outfile.write("%s, %2.1f, %s, %d, %s, " %
            #              (realization, MAGNITUDE, MECH, DIST, st_id))
            add_comma = 0
            for line in input_file:
                line = line.strip()
                # Skip comments
                if line.startswith("#"):
                    continue
                pieces = [float(item) for item in line.split()]
                # Write RotD50 value to file
                if add_comma:
                    outfile.write(",%10.7f" % (pieces[3]))
                else:
                    add_comma = 1
                    outfile.write("%10.7f" % (pieces[3]))
            # Done with this file
            input_file.close()
            # Add newline to outfile
            outfile.write("\n")
    # All done!
    outfile.close()

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

parser = optparse.OptionParser()
parser.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
parser.add_option("-o", "--output", dest="output_file",
                  help="Output plots' stem")
parser.add_option("-c", "--codebase", dest="codebase",
                  help="Codebase used in the simulation")
parser.add_option("-m", "--mag", type="float", dest="mag",
                  help="Event Magnitude")
parser.add_option("--dist", type="int", dest="dist",
                  help="Distance of the stations in km")
parser.add_option("--mech", dest="mech",
                  help="Fault mechanism")
parser.add_option("--vel", dest="vel",
                  help="Velocity model")
parser.add_option("--nga-west1", action="store_true", dest="nga1",
                  help="Uses NGA West 1 GMPEs for comparison")
(options, args) = parser.parse_args()


if options.input_dir is None:
    parser.error("Please specify the input directory!")
top_input_dir = options.input_dir
if not os.path.isdir(top_input_dir):
    parser.error("Invalid input directory!")
dirs = os.listdir(top_input_dir)
if not "Sims" in dirs:
    parser.error("Please provide the top-level simulation directory!\n"
                 "This is the directory given to the cluster script")
input_outdir = os.path.join(top_input_dir, "Sims", "outdata")
# Validate codebase to use
codebase = options.codebase
if codebase is None:
    print "Please specify a codebase!"
    sys.exit(1)
codebase = codebase.lower()
if codebase not in CODEBASES:
    print "Codebase needs to be one of: %s" % (CODEBASES)
    sys.exit(1)
if options.mag is None:
    parser.error("Please specify magnitude!")
MAGNITUDE = float(options.mag)
if options.mech is None:
    parser.error("Please specify fault mechanism!")
MECH = options.mech.upper()
if options.dist is None:
    parser.error("Please specify station distance!")
DIST = int(options.dist)
if options.vel is None:
    parser.error("Please specify velocity model!")
VEL = options.vel.upper()
if options.output_file is None:
    output_file = ("boxplot_%s-gmpe%d%s%d-%s" % (codebase,
                                                 int(MAGNITUDE * 10),
                                                 MECH,
                                                 DIST,
                                                 VEL))
else:
    output_file = options.output_file
# Default is to use NGA-WEST2
if options.nga1 is None:
    MODEL = "nga-west2"
else:
    MODEL = "nga-west1"

# Create temp dir
tmpdir = tempfile.mkdtemp(prefix="bbp-")
tmp_gmpe_file = os.path.join(tmpdir, "bbp_gmpe_data.txt")
tmp_sim_file = os.path.join(tmpdir, "bbp_sim_data.txt")

# Create data files with both gmpe and simulation data
create_gmpe_data_file(input_outdir, tmp_gmpe_file)
create_sim_data_file(input_outdir, tmp_sim_file)

# Run Matlab script
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
bband_utils.runprog('export MATLABPATH=%s; ' % (SCRIPT_PATH) +
                    'matlab -nodisplay -r ' +
                    '"SimsFile=\'%s\';' % (tmp_sim_file) +
                    'GMPEFile=\'%s\';' % (tmp_gmpe_file) +
                    'OUTFile=\'%s\';' % (output_file) +
                    'Method=\'%s\';Mag=%f;' % (codebase, MAGNITUDE) +
                    'Vel=\'%s\';' % (VEL) +
                    'Model=\'%s\';' % (MODEL) +
                    'Dist=%d;Mech=\'%s\';%s"' % (DIST, MECH, BOXPLOT) +
                    '>/dev/null; stty echo')
print "All Done!"
# Clean-up, all done!
shutil.rmtree(tmpdir)
