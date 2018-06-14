#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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

    args.general_record_seq_no = "-999"
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

def write_output_data2(args):
    """
    This function writes all output to the flat file
    """
    # Output filenames
    output_filename_prefix = args.prefix
    output_filename_suffix = args.suffix
    if args.suffix:
        output_filename_suffix = "_%s" % (args.suffix)
    output_filename_date = datetime.date.today().strftime("%y%m%d")
    output_filename_extension = ".f01.csv"
    output_main_filename = os.path.join(args.output_dir,
                                        "%s_%s%s_Main%s" %
                                        (output_filename_prefix,
                                         output_filename_date,
                                         output_filename_suffix,
                                         output_filename_extension))
    output_psa_h1_filename = os.path.join(args.output_dir,
                                          "%s_%s%s_PSA_H1_D0pt05%s" %
                                          (output_filename_prefix,
                                           output_filename_date,
                                           output_filename_suffix,
                                           output_filename_extension))
    output_psa_h2_filename = os.path.join(args.output_dir,
                                          "%s_%s%s_PSA_H2_D0pt05%s" %
                                          (output_filename_prefix,
                                           output_filename_date,
                                           output_filename_suffix,
                                           output_filename_extension))
    output_psa_v_filename = os.path.join(args.output_dir,
                                         "%s_%s%s_PSA_V_D0pt05%s" %
                                         (output_filename_prefix,
                                          output_filename_date,
                                          output_filename_suffix,
                                          output_filename_extension))
    output_psa_rd50_filename = os.path.join(args.output_dir,
                                            "%s_%s%s_PSA_RotD50_D0pt05%s" %
                                            (output_filename_prefix,
                                             output_filename_date,
                                             output_filename_suffix,
                                             output_filename_extension))
    output_psa_rd100_filename = os.path.join(args.output_dir,
                                             "%s_%s%s_PSA_RotD100_D0pt05%s" %
                                             (output_filename_prefix,
                                              output_filename_date,
                                              output_filename_suffix,
                                              output_filename_extension))

    # Create header
    header = ("bbp_software_version,sim_simulation_workflow,"
              "sim_method_short_name,sim_site_effects,"
              "record_sequence_number,eq_id,eq_magnitude,"
              "realization,number_of_segments")
    header = ("%s,segment_lengths,segment_widths,segment_ztors,"
              "segment_strikes,segment_rakes,segment_dips,"
              "total_length,average_strike,average_rake,"
              "average_dip,average_width,average_ztor,"
              "mechanism_based_on_average_rake" % (header))
    header = ("%s,vmodel_name,gf_name,gf_dt,vmodel_vs30" % (header))
    header = ("%s,sim_station_name,sim_station_latitude,"
              "sim_station_longitude,sim_station_elevation,"
              "target_station_vs30,"
              "target_station_nehrp_class,station_rrup,station_rjb,station_rx,"
              "num_components,acc_file_name,"
              "h1_azimuth,h2_azimuth,v_orientation,"
              "dt,num_samples,nyquist,luf,huf,ufb,lup,hup,upb,"
              "pga_h1,pga_h2,pga_v,pgv_h1,pgv_h2,pgv_v" % (header))
    header = ("%s,ai_h1,ai_h2,ai_v,aid5_75_h1,aid5_75_h2,aid5_75_v,"
              "aid5_95_h1,aid5_95_h2,aid5_95_v" % (header))

    header_psa = "record_sequence_number,intensity_measure,damping"
    for period in args.rd50_periods:
        header_psa = ("%s,T%dp%03d" % (header_psa,
                                       int(period),
                                       (period % 1 * 1000)))

    # Create first (common) part of the output
    sim_params = ('"%s","%s","%s","%s","%s","%s",%s' %
                  (args.bbp_software_info_version,
                   "/".join(args.bbp_software_info_modules),
                   args.general_method,
                   args.bbp_software_info_site,
                   args.general_record_seq_no,
                   args.general_eqid,
                   str(args.general_magnitude)))

    # Output main data file
    output_file = open(output_main_filename, 'w')
    # Print header
    output_file.write('%s\n' % (header))
    for realization in args.realizations:
        realization_data = args.data[realization]
        station_names = realization_data["station_names"]
        realization_params = ('%s,%d' % (realization,
                                         realization_data["num_src"]))
        realization_params = ('%s,"%s","%s","%s","%s","%s","%s",%s,%s,%s,'
                              '%s,%s,%s,"%s"' %
                              (realization_params,
                               realization_data["segments_length"],
                               realization_data["segments_width"],
                               realization_data["segments_ztor"],
                               realization_data["segments_strike"],
                               realization_data["segments_rake"],
                               realization_data["segments_dip"],
                               realization_data["total_length"],
                               realization_data["average_strike"],
                               realization_data["average_rake"],
                               realization_data["average_dip"],
                               realization_data["average_width"],
                               realization_data["average_ztor"],
                               realization_data["average_mechanism"]))
        realization_params = ('%s,"%s","%s",%.2f,%s' %
                              (realization_params,
                               realization_data["vmodel_name"],
                               realization_data["gf_name"],
                               realization_data["gf_dt"],
                               realization_data["vs_30"]))
        for station in station_names:
            st_data = realization_data["stations"][station]
            station_params = ('%s,%s,%s,%.1f,%s,"%s",%s,%s,%s,%s,'
                              '"%s",%s,%s,"%s",%s,%s,%s,%s,%s,%s,'
                              '%s,%s,%s,%s,%s,%s,%s,%s,%s' %
                              (station,
                               st_data["sim_station_latitude"],
                               st_data["sim_station_longitude"],
                               st_data["sim_station_elevation"],
                               st_data["target_station_vs30"],
                               st_data["target_station_nehrp_class"],
                               st_data["rrup"],
                               st_data["rjb"],
                               st_data["rx"],
                               st_data["components"],
                               st_data["acc_file_name"],
                               st_data["h1_azimuth"],
                               st_data["h2_azimuth"],
                               st_data["v_orientation"],
                               st_data["time_series_dt"],
                               st_data["time_series_num_samples"],
                               st_data["nyquist"],
                               st_data["luf"],
                               st_data["huf"],
                               st_data["ufb"],
                               st_data["lup"],
                               st_data["hup"],
                               st_data["upb"],
                               st_data["pga_h1"],
                               st_data["pga_h2"],
                               st_data["pga_v"],
                               st_data["pgv_h1"],
                               st_data["pgv_h2"],
                               st_data["pgv_v"]))
            station_params = ('%s,%.2f,%.2f,%.2f,'
                              '%.2f,%.2f,%.2f,'
                              '%.2f,%.2f,%.2f' % (station_params,
                                                  st_data["ai_h1"],
                                                  st_data["ai_h2"],
                                                  st_data["ai_v"],
                                                  st_data["ad5_75_h1"],
                                                  st_data["ad5_75_h2"],
                                                  st_data["ad5_75_v"],
                                                  st_data["ad5_95_h1"],
                                                  st_data["ad5_95_h2"],
                                                  st_data["ad5_95_v"]))
            output_file.write('%s,%s,%s\n' %
                              (sim_params, realization_params, station_params))

    # All done
    output_file.close()

    # Write PSA files
    psa_files = [output_psa_h1_filename, output_psa_h2_filename,
                 output_psa_v_filename, output_psa_rd50_filename,
                 output_psa_rd100_filename]
    psa_measurements = ["psa_h1", "psa_h2", "psa_v",
                        "rd50", "rd100"]
    psa_meas_labels = ["PSA_H1", "PSA_H2", "PSA_V",
                       "PSA_RotD50", "PSA_RotD100"]

    for output_filename, psa_data, psa_label in zip(psa_files,
                                                    psa_measurements,
                                                    psa_meas_labels):
        # Output psa data file
        output_file = open(output_filename, 'w')
        # Print header
        output_file.write("%s\n" % (header_psa))

        for realization in args.realizations:
            realization_data = args.data[realization]
            station_names = realization_data["station_names"]
            for station in station_names:
                st_data = realization_data["stations"][station]
                if st_data[psa_data] is None:
                    continue
                psa_params = '"%s","%s",%.2f' % (st_data["acc_file_name"],
                                                 psa_label,
                                                 st_data["damping"])
                for period in st_data[psa_data]:
                    psa_params = ('%s,%.7f' % (psa_params, period))

                # Write output
                output_file.write('%s\n' % (psa_params))

        # All done
        output_file.close()

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
