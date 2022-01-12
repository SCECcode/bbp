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

This module takes care of building a workflow using either user
choices interactively, or an option file containing all needed
parameters.

This module compiles results from a single simulation and creates
a csv file containing information about all stations.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import argparse
import datetime
import numpy as np
from scipy import integrate

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList
import bband_utils
import xml_handler
import velocity_models
from export_bbp_cluster_simulation import calculate_vs30, \
    calculate_mechanism, read_bbp, calculate_arias, \
    calculate_distances, get_vmodel_from_html, get_method_from_html, \
    collect_rd50_values, collect_rd100_values, \
    collect_station_params, write_output_data

# Initialize global variables
INSTALL = InstallCfg.getInstance()
CM2G = 980.664999

def parse_arguments():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Exports data from a BBP "
                                     "simulation to a flat file.")
    parser.add_argument("--sim_id", dest="sim_id",
                        type=int, required=True,
                        help="simulation id")
    parser.add_argument("--output_dir", "-o", dest="output_dir",
                        required=True, help="output directory")
    parser.add_argument("--copy", "-c", dest="copy_timeseries", action='store_true',
                        help="copy all timeseries into the output directory")
    parser.add_argument("--prefix", "-p", dest="prefix",
                        default="BBP_Study",
                        help="prefix for output files")
    parser.add_argument("--suffix", "-s", dest="suffix",
                        default="",
                        help="suffix for output files")
    parser.add_argument("--first_seg_id", dest="first_seg_id", type=int,
                        help="sim_id of first segment for multisegment simulations")
    args = parser.parse_args()

    # Default value is sim_id
    if args.first_seg_id is None:
        args.first_seg_id = args.sim_id

    return args

def collect_simulation_params(args):
    """
    This function collects simulation-wide parameters
    """
    # Get paths to one xml and a SRC file
    xml_path = os.path.join(INSTALL.A_XML_DIR,
                            "%s.xml" % (str(args.first_seg_id)))
    src_files = glob.glob("%s/*.src" % (args.top_level_indir))
    src_path = os.path.join(args.top_level_indir, src_files[0])
    html_file = glob.glob("%s/*.html" % (args.top_level_htmldir))[0]

    # Get simulation method from html file
    args.general_method = get_method_from_html(html_file).lower()

    # Parse SRC and get magnitude
    src_keys = bband_utils.parse_src_file(src_path)
    args.general_magnitude = src_keys["magnitude"]

    # Parse XML file
    workflow_obj = xml_handler.parse_xml(xml_path)
    args.bbp_software_info_version = str(workflow_obj.version)
    modules = []
    for item in workflow_obj.workflow:
        modules.append(str(item.getName()))
    args.bbp_software_info_modules = modules
    if "WccSiteamp" in modules:
        args.bbp_software_info_site = "GP2014"
    else:
        args.bbp_software_info_site = "None"

    args.general_eqid = "-999"

def collect_realization_params(args, realization):
    """
    Collects parameters for one realization
    """
    src_files = glob.glob("%s/*.src" % (args.top_level_indir))
    stl_file = glob.glob("%s/*.stl" % (args.top_level_indir))[0]
    data = {}

    # Compile data from SRC file(s)
    data["num_src"] = len(src_files)
    # Save info in args too for first realization
    if "num_src" not in args:
        args.num_src = len(src_files)
    for i, src_file in zip(range(1, len(src_files) + 1), src_files):
        src_index = "bbp_src_%d" % (i)
        src_keys = bband_utils.parse_src_file(src_file)
        src_keys["mechanism"] = calculate_mechanism(src_keys["rake"])
        data[src_index] = src_keys

    # Combine SRC information
    data["segments_length"] = data["bbp_src_1"]["fault_length"]
    data["segments_width"] = data["bbp_src_1"]["fault_width"]
    data["segments_ztor"] = data["bbp_src_1"]["depth_to_top"]
    data["segments_strike"] = data["bbp_src_1"]["strike"]
    data["segments_rake"] = data["bbp_src_1"]["rake"]
    data["segments_dip"] = data["bbp_src_1"]["dip"]
    data["total_length"] = float(data["bbp_src_1"]["fault_length"])
    data["average_strike"] = [float(data["bbp_src_1"]["strike"])]
    data["average_rake"] = [float(data["bbp_src_1"]["rake"])]
    data["average_dip"] = [float(data["bbp_src_1"]["dip"])]
    data["average_width"] = [float(data["bbp_src_1"]["fault_width"])]
    data["average_ztor"] = [float(data["bbp_src_1"]["depth_to_top"])]
    for i in range(2, len(src_files) + 1):
        src_index = "bbp_src_%d" % (i)
        data["segments_length"] = "%s,%s" % (data["segments_length"],
                                             data[src_index]["fault_length"])
        data["segments_width"] = "%s,%s" % (data["segments_width"],
                                            data[src_index]["fault_width"])
        data["segments_ztor"] = "%s,%s" % (data["segments_ztor"],
                                           data[src_index]["depth_to_top"])
        data["segments_strike"] = "%s,%s" % (data["segments_strike"],
                                             data[src_index]["strike"])
        data["segments_rake"] = "%s,%s" % (data["segments_rake"],
                                             data[src_index]["rake"])
        data["segments_dip"] = "%s,%s" % (data["segments_dip"],
                                          data[src_index]["dip"])
        data["total_length"] = (data["total_length"] +
                                float(data[src_index]["fault_length"]))
        data["average_strike"].append(data[src_index]["strike"])
        data["average_rake"].append(data[src_index]["rake"])
        data["average_dip"].append(data[src_index]["dip"])
        data["average_width"].append(data[src_index]["fault_width"])
        data["average_ztor"].append(data[src_index]["depth_to_top"])
    data["average_strike"] = np.average(data["average_strike"])
    data["average_rake"] = np.average(data["average_rake"])
    data["average_dip"] = np.average(data["average_dip"])
    data["average_width"] = np.average(data["average_width"])
    data["average_ztor"] = np.average(data["average_ztor"])
    data["average_mechanism"] = calculate_mechanism(data["average_rake"])

    # Get velocity model data
    html_file = glob.glob("%s/*.html" % (args.top_level_htmldir))[0]
    data["vmodel_name"] = get_vmodel_from_html(html_file)
    vel_obj = velocity_models.get_velocity_model_by_name(data["vmodel_name"])
    if vel_obj is None:
        print("ERROR: Cannot find velocity model %s!" % (data["vmodel_name"]))
        sys.exit(-1)
    if args.general_method in ["gp", "sdsu", "song"]:
        vmodel_params = vel_obj.get_codebase_params('gp')
        vmodel_file = vel_obj.get_velocity_model('gp')
        data["gf_name"] = vmodel_params['GF_NAME']
        data["vs_30"] = calculate_vs30(vmodel_file)
        data["gf_dt"] = float(vmodel_params['GF_DT'])
    elif args.general_method in ["ucsb"]:
        vmodel_params = vel_obj.get_codebase_params('ucsb')
        vmodel_file = vel_obj.get_velocity_model('ucsb')
        data["gf_name"] = vmodel_params['GREEN_SOIL']
        data["vs_30"] = "-999"
        data["gf_dt"] = float(vmodel_params['GF_DT'])
    else:
        data["gf_name"] = "-888"
        data["vs_30"] = "-888"
        data["gf_dt"] = "-888"

    # Parse STL file
    slo = StationList(stl_file)
    site_list = slo.getStationList()
    station_names = []
    for site in site_list:
        station_names.append(site.scode)
    data["station_names"] = station_names

    stations = {}
    for site in site_list:
        stations[site.scode] = {}
        if args.bbp_software_info_site == "None":
            vs_30 = data["vs_30"]
        elif site.vs30 is None:
            vs_30 = data["vs_30"]
        else:
            vs_30 = site.vs30
        collect_station_params(site, stations[site.scode], src_files,
                               args, realization, vs_30)
        collect_rd50_values(stations[site.scode], args)
        collect_rd100_values(stations[site.scode], args)

    # Save data
    data["stations"] = stations

    # Save realization data
    args.data[realization] = data

def create_flat_file_from_bbp_sim():
    """
    Create a flat file from a BBP simulation
    """
    # Get all we need from the command-line
    args = parse_arguments()

    # Figure out top-level directories
    args.top_level_indir = os.path.join(INSTALL.A_IN_DATA_DIR,
                                        str(args.sim_id))
    args.top_level_outdir = INSTALL.A_OUT_DATA_DIR
    args.top_level_htmldir = os.path.join(INSTALL.A_OUT_DATA_DIR,
                                          str(args.first_seg_id))
    args.realizations = [str(args.sim_id)]
    args.data = {}

    # Create top-level output directory
    bband_utils.mkdirs([args.output_dir], print_cmd=False)

    # Collect simulation-wide parameters
    collect_simulation_params(args)

    # Collect parameters for each realization
    print("==> Processing simulation : %d..." % (args.sim_id))
    # Create output directory for realization if requested to copy seismograms
    if args.copy_timeseries:
        bband_utils.mkdirs([os.path.join(args.output_dir,
                                         str(args.sim_id))],
                           print_cmd=False)
    collect_realization_params(args, str(args.sim_id))

    # Write flat file
    write_output_data(args)

if __name__ == '__main__':
    create_flat_file_from_bbp_sim()
