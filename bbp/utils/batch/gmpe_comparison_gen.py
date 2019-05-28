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

Python program to generate PSAs using the various NGA methods for any
simulation. This code uses Pynga to calculate the PSAs using the SRC
file and a station list. It then compares the generated results
against PSAs from synthetic seismograms.
"""

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
import gmpe_config
from station_list import StationList

# Import Pynga and its utilities
import pynga
import pynga.utils as putils

# --------------------------------------------------------------------------
# Some constants
# --------------------------------------------------------------------------
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "irikura"]
PERIODS = [0.010, 0.011, 0.012, 0.013, 0.015, 0.017, 0.020, 0.022, 0.025,
           0.029, 0.032, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065,
           0.075, 0.085, 0.100, 0.110, 0.120, 0.130, 0.150, 0.170, 0.200,
           0.220, 0.240, 0.260, 0.280, 0.300, 0.350, 0.400, 0.450, 0.500,
           0.550, 0.600, 0.650, 0.750, 0.850, 1.000, 1.100, 1.200, 1.300,
           1.500, 1.700, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, 3.500,
           4.000, 4.400, 5.000, 5.500, 6.000, 6.500, 7.500, 8.500, 10.000]
BOXPLOT = "gmpe_comparison"

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def parse_src_file(a_srcfile):
    """
    Function parses the SRC file and checks for needed keys. It
    returns a dictionary containing the keys found in the src file.
    """
    src_keys = bband_utils.parse_properties(a_srcfile)
    required_keys = ["magnitude", "fault_length", "fault_width", "dlen",
                     "dwid", "depth_to_top", "strike", "rake", "dip",
                     "lat_top_center", "lon_top_center"]
    for key in required_keys:
        if key not in src_keys:
            raise bband_utils.ParameterError("key %s missing in src file" %
                                             (key))
    # Convert keys to floats
    for key in src_keys:
        src_keys[key] = float(src_keys[key])

    return src_keys

def calculate_gmpe(src_keys, station, output_file, rrups, gmpe_group_name):
    """
    Calculate the GMPE results for a given station.
    """
    gmpe_group = gmpe_config.GMPES[gmpe_group_name]
    origin = (src_keys['lon_top_center'], src_keys['lat_top_center'])
    dims = (src_keys['fault_length'], src_keys['dlen'],
            src_keys['fault_width'], src_keys['dwid'],
            src_keys['depth_to_top'])
    mech = (src_keys['strike'], src_keys['dip'], src_keys['rake'])

    # Station location
    site_geom = [float(station.lon), float(station.lat), 0.0]
    (fault_trace1, upper_seis_depth,
     lower_seis_depth, ave_dip,
     dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
    rjb, rrup, rx = putils.DistanceToSimpleFaultSurface(site_geom,
                                                        fault_trace1,
                                                        upper_seis_depth,
                                                        lower_seis_depth,
                                                        ave_dip)

    print "station: %s, Rrup: %f" % (station.scode, rrup)
    rrups.append(rrup)

    vs30 = 1000
    z10 = None # Let PyNGA calculate it
    z25 = None # Let PyNGA calculate it

    # Compute PSA for this stations
    station_median = []
    for period in gmpe_group["periods"]:
        period_medians = []
        for nga_model in gmpe_group["models"]:
            median = gmpe_config.calculate_gmpe(gmpe_group_name,
                                                nga_model,
                                                src_keys['magnitude'],
                                                rjb, vs30,
                                                period,
                                                rake=src_keys['rake'],
                                                dip=src_keys['dip'],
                                                W=src_keys['fault_width'],
                                                Ztor=src_keys['depth_to_top'],
                                                Rrup=rrup, Rx=rx,
                                                Z10=z10, Z25=z25)
            period_medians.append(median)
        station_median.append((period, period_medians))

    # Create label
    file_label = ""
    for nga_model in gmpe_group["models"]:
        file_label = "%s %s" % (file_label, nga_model)
    # Output data to file
    outfile = open(output_file, 'w')
    outfile.write("#station: %s\n" % (station.scode))
    outfile.write("#period%s\n" % (file_label))
    for item in station_median:
        period = item[0]
        vals = item[1]
        out_str = "%.4f" % (period)
        for method in vals:
            out_str = out_str + "\t%.6f" % (method)
        outfile.write("%s\n" % (out_str))
    outfile.close()

    # Return list
    return station_median

def create_gmpe_data_file(indata_dir, tmpdir,
                          gmpe_file, gmpe_label_file,
                          gmpe_group_name):
    """
    This function creates a file containing the GMPE data
    """
    # Find SRC file
    basedir = os.path.join(indata_dir, os.listdir(indata_dir)[0])
    src_file = glob.glob("%s%s*.src" % (basedir, os.sep))
    if not len(src_file):
        print "Unable to find SRC file!"
        sys.exit(1)
    src_file = src_file[0]
    # Now parse SRC file
    src_keys = parse_src_file(src_file)

    # Find station list
    stl_file = glob.glob("%s%s*.stl" % (basedir, os.sep))
    if len(stl_file) != 1:
        print "Unable to find STL file!"
        sys.exit(1)
    stl_file = stl_file[0]
    # Parse station list
    slo = StationList(stl_file)
    site_list = slo.getStationList()

    # Write ri50 files
    rrups = []
    for site in site_list:
        output_file = os.path.join(tmpdir, "%s.ri50" % (site.scode))
        calculate_gmpe(src_keys, site, output_file, rrups, gmpe_group_name)
    mean_rrup = numpy.mean(rrups)

    # Get periods
    gmpe_group = gmpe_config.GMPES[gmpe_group_name]

    # Write label file
    out_labels = open(gmpe_label_file, 'w')
    # Write labels
    labels = ",".join(gmpe_group["labels"])
    out_labels.write("%s\n" % (labels))
    # Done
    out_labels.close()

    # Open output file, write header
    outfile = open(gmpe_file, 'w')
    # Add header for the GMPE column
    outfile.write("0")
    for period in gmpe_group["periods"]:
        outfile.write(",%10.5f" % period)
    outfile.write("\n")

    # Get number of GMPEs that we have
    number_of_gmpes = len(gmpe_group["models"])

    # Get list of stations to process
    stations = sorted(glob.glob("%s%s*.ri50" % (tmpdir, os.sep)))
    for station in stations:
        # Start empty
        gmpe_ri50 = []

        input_file = open(station, 'r')
        for line in input_file:
            line = line.strip()
            # Skip comments
            if line.startswith("#"):
                continue
            pieces = [float(item) for item in line.split()]
            # Initialize gmpe_ri50 structure
            if not gmpe_ri50:
                for item in pieces[1:]:
                    gmpe_ri50.append([])
            for item, dst in zip(pieces[1:], gmpe_ri50):
                dst.append(item)
        # Done with input file
        input_file.close()
        # Read all values
        for i in range(0, len(gmpe_ri50)):
            outfile.write("%d" % (i + 1))
            for item in gmpe_ri50[i]:
                outfile.write(",%10.6f" % (item))
            outfile.write("\n")

    # All done, close output file
    outfile.close()

    return (src_keys['magnitude'], mean_rrup, number_of_gmpes)

def create_sim_data_file(outdata_dir, sim_file):
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

    realizations = sorted(os.listdir(outdata_dir))
    for realization in realizations:
        basedir = os.path.join(outdata_dir, realization)
        stations = sorted(glob.glob("%s%s*.rd50" % (basedir, os.sep)))
        for station in stations:
            #st_id = os.path.basename(station)
            # Find just the numeric part of the station id
            #st_id = st_id[st_id.find('-')+1:st_id.find(".rd50")]
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

gmpe_groups_available = gmpe_config.GMPES.keys()
gmpe_groups_available_lc = [gmpe.lower() for gmpe in gmpe_groups_available]

parser = optparse.OptionParser()
parser.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
parser.add_option("-o", "--output", dest="output_file",
                  help="Output plots' stem")
parser.add_option("-c", "--codebase", dest="codebase",
                  help="Codebase used in the simulation")
parser.add_option("-l", "--label", dest="label",
                  help="Label for the scenario or validation simulation")
parser.add_option("-g", "--gmpe-group", dest="gmpe_group_name",
                  help="GMPE group: %s" % (gmpe_groups_available_lc))
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
outdata_dir = os.path.join(top_input_dir, "Sims", "outdata")
indata_dir = os.path.join(top_input_dir, "Sims", "indata")
# Validate codebase to use
codebase = options.codebase
if codebase is None:
    print "Please specify a codebase!"
    sys.exit(1)
codebase = codebase.lower()
if codebase not in CODEBASES:
    print "Codebase needs to be one of: %s" % (CODEBASES)
    sys.exit(1)
label = options.label
if label is None:
    print "Please specify simulation label!"
    sys.exit(1)
if options.output_file is None:
    output_file = ("gmpe-comparison-%s-%s" % (codebase, label))
else:
    output_file = options.output_file

# Figure out which gmpe group to use
gmpe_group_name = options.gmpe_group_name
if gmpe_group_name is None:
    print "Please specify gmpe group name!"
    sys.exit(1)
if not gmpe_group_name.lower() in gmpe_groups_available_lc:
    print "Invalid gmpe group name!"
    print "Options are: %s" % (gmpe_groups_available_lc)
    sys.exit(1)
gmpe_group_index = gmpe_groups_available_lc.index(gmpe_group_name.lower())
gmpe_group_name = gmpe_groups_available[gmpe_group_index]

# Create temp dir
tmpdir = tempfile.mkdtemp(prefix="bbp-")
tmp_gmpe_file = os.path.join(tmpdir, "bbp_gmpe_data.txt")
tmp_gmpe_label_file = os.path.join(tmpdir, "bbp_gmpe_labels.txt")
tmp_sim_file = os.path.join(tmpdir, "bbp_sim_data.txt")

# Create data files with both gmpe and simulation data
(mag, mean_rrup, num_gmpes) = create_gmpe_data_file(indata_dir, tmpdir,
                                                    tmp_gmpe_file,
                                                    tmp_gmpe_label_file,
                                                    gmpe_group_name)
create_sim_data_file(outdata_dir, tmp_sim_file)

# Run Matlab script
SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
bband_utils.runprog('export MATLABPATH=%s; ' % (SCRIPT_PATH) +
                    'matlab -nodisplay -r ' +
                    '"SimsFile=\'%s\';' % (tmp_sim_file) +
                    'GMPEFile=\'%s\';' % (tmp_gmpe_file) +
                    'GMPELabels=\'%s\';' % (tmp_gmpe_label_file) +
                    'OUTFile=\'%s\';' % (output_file) +
                    'NumGMPE=%d;' % (num_gmpes) +
                    'PlotTitle=\'%s, %s, M=%.2f, Mean(Rrup)=%.2f\';' %
                    (label, codebase.upper(), mag, mean_rrup) +
                    '%s"' % (BOXPLOT) +
                    '>/dev/null; stty echo')
print "All Done!"
# Clean-up, all done!
shutil.rmtree(tmpdir)
