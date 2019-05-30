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

Program to plot a GoF comparison between two sets of seismograms
"""

# Import Python modules
import os
import sys
import glob
import shutil
import tempfile

# Import Broadband modules
import bband_utils
from gp_gof_cfg import GPGofCfg
from install_cfg import InstallCfg
from station_list import StationList
from PlotGOF import PlotGoF

# Import Pynga and its utilities
import pynga.utils as putils

# Main
if len(sys.argv) != 6:
    print ("Usage: %s station_list src_file sim_id_1 sim_id_2 output_dir" %
           (os.path.basename(sys.argv[0])))
    sys.exit(1)

# Create temp dir
TMPDIR = tempfile.mkdtemp(prefix="bbp-")
resid_file = os.path.join(TMPDIR, "bbp-rd50-resid.txt")
log_file = os.path.join(TMPDIR, "bbp-rd50-resid.log")

# Get input parameters
station_list = sys.argv[1]
src_file = sys.argv[2]
sim_id_1 = int(sys.argv[3])
sim_id_2 = int(sys.argv[4])
output_dir = sys.argv[5]

# Create directory paths
install = InstallCfg.getInstance()
config = GPGofCfg()
a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id_1))
a_outdir1 = os.path.join(install.A_OUT_DATA_DIR, str(sim_id_1))
a_outdir2 = os.path.join(install.A_OUT_DATA_DIR, str(sim_id_2))

# Src file
a_srcfile = os.path.join(a_indir, src_file)
src_keys = bband_utils.parse_src_file(a_srcfile)

# Station file
a_statfile = os.path.join(a_indir, station_list)
slo = StationList(a_statfile)
site_list = slo.getStationList()

# Capture event_label
bias_file = glob.glob("%s%s*.bias" % (a_outdir1, os.sep))
if len(bias_file) < 1:
    raise bband_utils.ProcessingError("Cannot find event label!")
bias_file = bias_file[0]
# Let's capture the event label
event_label = os.path.basename(bias_file).split("-")[0]

print_header_rd50 = 1

# Go through the stations
for site in site_list:
    stat = site.scode
    slon = float(site.lon)
    slat = float(site.lat)

    # Check files are there
    a_rd50file1 = os.path.join(a_outdir1, "%d.%s.rd50" %
                               (sim_id_1, stat))
    a_rd50file2 = os.path.join(a_outdir2, "%d.%s.rd50" %
                               (sim_id_2, stat))
    if not os.path.exists(a_rd50file1) or not os.path.exists(a_rd50file2):
        # Just skip it
        print "Skipping station %s..." % (stat)
        continue

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

    bband_utils.check_path_lengths([a_rd50file1, a_rd50file2, resid_file],
                                   bband_utils.GP_MAX_FILENAME)
    cmd = ("%s/gen_resid_tbl_3comp bbp_format=1 " %
           (install.A_GP_BIN_DIR) +
           "datafile1=%s simfile1=%s " % (a_rd50file1, a_rd50file2) +
           "comp1=psa5n comp2=psa5e comp3=rotd50 " +
           "eqname=%s mag=%s stat=%s lon=%.4f lat=%.4f " %
           (event_label, src_keys['magnitude'], stat, slon, slat) +
           "vs30=%d cd=%.2f " % (site.vs30, rrup) +
           "flo=%f fhi=%f " % (site.low_freq_corner,
                               site.high_freq_corner) +
           "print_header=%d >> %s 2>> %s" %
           (print_header_rd50, resid_file, log_file))
    bband_utils.runprog(cmd, abort_on_error=True)

    # Only need to print header the first time
    if print_header_rd50 == 1:
        print_header_rd50 = 0

for comp in config.COMPS_PSA5:
    # Build paths and check lengths
    fileroot = os.path.join(TMPDIR, "%s-%d_r%d-all-rd50-%s" %
                            (event_label, sim_id_1, config.MIN_CDST,
                             comp))
    bband_utils.check_path_lengths([fileroot],
                                   bband_utils.GP_MAX_FILENAME)

    cmd = ("%s/resid2uncer_varN " % (install.A_GP_BIN_DIR) +
           "residfile=%s fileroot=%s " % (resid_file, fileroot) +
           "comp=%s nstat=%d nper=63 " % (comp, len(site_list)) +
           "min_cdst=%d >> %s 2>&1" %
           (config.MIN_CDST, log_file))
    bband_utils.runprog(cmd, abort_on_error=True)

plottitle = ("GOF Comparison between simulation %d and simulation %d" %
             (sim_id_1, sim_id_2))
plot_mode = 'rd50-single'
fileroot = '%s-%d_r0-all-rd50' % (event_label, sim_id_1)
plotter = PlotGoF()
plotter.plot(plottitle, fileroot, TMPDIR, output_dir,
             0, plot_mode, 'single')

print "All Done!"
# Clean-up, all done!
shutil.rmtree(TMPDIR)
