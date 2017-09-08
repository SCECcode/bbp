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

This module takes care of building a workflow using either user
choices interactively, or an option file containing all needed
parameters.

This module compiles results from a cluster simulation and creates
a csv file containing information about all stations/realizations.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import shutil
import argparse

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList
import bband_utils
import xml_handler
import velocity_models

# Import Pynga and its utilities
import pynga.utils as putils

# Initialize global variables
INSTALL = InstallCfg.getInstance()

def parse_arguments():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Exports data from a BBP cluster "
                                     "simulation to a flat file.")
    parser.add_argument("--input_dir", "-i", dest="input_dir",
                        help="input directory")
    parser.add_argument("--output_dir", "-o", dest="output_dir",
                        help="output directory")
    parser.add_argument("--copy", "-c", dest="copy_timeseries", action='store_true',
                        help="copy all timeseries into the output directory")
    args = parser.parse_args()

    if args.input_dir is None:
        print("[ERROR]: Please provide input directory!")
        sys.exit(-1)
    else:
        if not os.path.isdir(args.input_dir):
            print("[ERROR]: Please provide valid input directory!")
            sys.exit(-1)
        if not "Sims" in os.listdir(args.input_dir):
            print("[ERROR]: Please provide top-level cluster simulation directory!")
            sys.exit(-1)
    if args.output_dir is None:
        print("[ERROR]: Please provide output directory!")
        sys.exit(-1)

    return args

def collect_simulation_params(args):
    """
    This function collects simulation-wide parameters
    """
    # Get paths to one xml and a SRC file
    first_realization = args.realizations[0]
    xml_dir = os.path.join(args.input_dir, "Xml")
    xml_files = glob.glob("%s/*.xml" % (xml_dir))
    xml_path = os.path.join(xml_dir, xml_files[0])
    src_dir = os.path.join(args.top_level_indir, args.realizations[0])
    src_files = glob.glob("%s/*.src" % (src_dir))
    src_path = os.path.join(src_dir, src_files[0])
    html_dir = os.path.join(args.top_level_outdir, args.realizations[0])
    html_file = glob.glob("%s/*.html" % (html_dir))[0]

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
        
    if workflow_obj.val_obj is not None:
        val_obj = workflow_obj.val_obj
        args.general_eq_name = workflow_obj.val_obj.get_print_name()
    else:
        args.general_eq_name = "NA"

    args.general_earthquake_location = "TBD"
    args.general_record_seq_no = "TBD"
    args.general_eqid = "TBD"
    args.general_fault_information_source = "TBD"
    args.general_fault_id = "TBD"
    args.general_fault_name = "TBD"

def calculate_mechanism(rake):
    """
    Compute Mechanism based on rake angle
    """
    if rake >= -180 and rake < -150:
        return "Strike-Slip"
    if rake >= -30 and rake <= 30:
        return "Strike-Slip"
    if rake > 150 and rake <= 180:
        return "Strike-Slip"
    if rake >= -120 and rake < -60:
        return "Normal"
    if rake > 60 and rake <= 120:
        return "Reverse"
    if rake > 30 and rake <= 60:
        return "Reverse-Oblique"
    if rake > 120 and rake <= 150:
        return "Reverse-Oblique"
    if rake >= -150 and rake < -120:
        return "Normal-Oblique"
    if rake >= -60 and rake < -30:
        return "Normal-Oblique"
    return "Unknown"

def get_vmodel_from_html(html_file):
    """
    Parse vmodel name from html_file
    """
    input_file = open(html_file, 'r')
    for line in input_file:
        if line.find("Velocity model version") > 0:
            line = next(input_file)
            break
    token = line[4:].split(" ")[0]
    input_file.close()

    return token

def get_method_from_html(html_file):
    """
    Parse simulation method name from html_file
    """
    input_file = open(html_file, 'r')
    for line in input_file:
        if line.find("Simulation Method") > 0:
            line = next(input_file)
            break
    token = line[4:].split("<")[0]
    input_file.close()

    return token

def calculate_distances(src_files, site):
    """
    Calculate Rrup, Rjb, Rx using multiple SRC files
    """
    rrup = 10000000
    rjb = 10000000
    rx = 10000000
    
    for src_file in src_files:
        src_keys = bband_utils.parse_src_file(src_file)        
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
        my_rjb, my_rrup, my_rx = putils.DistanceToSimpleFaultSurface(site_geom,
                                                                     fault_trace1,
                                                                     up_seis_depth,
                                                                     low_seis_depth,
                                                                     ave_dip)
        rjb = min(my_rjb, rjb)
        rrup = min(my_rrup, rrup)
        rx = min(my_rx, rx)        

    return rrup, rjb, rx

def calculate_timeseries_param(station, site, args, realization):
    """
    Calculate/collect parameters from timeseries
    """
    vel_file = os.path.join(args.top_level_outdir, station["vel_file_name"])
    acc_file = os.path.join(args.top_level_outdir, station["acc_file_name"])

    # Read velocity timeseries
    num_samples = 0
    time_0 = -999
    time_1 = -999
    dt = -999
    pgv1 = -999
    pgv2 = -999
    pgv3 = -999
    pga1 = -999
    pga2 = -999
    pga3 = -999
    
    input_bbp = open(vel_file)
    for line in input_bbp:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        num_samples = num_samples + 1
        if dt == -999:
            if time_0 == -999:
                time_0 = tokens[0]
            elif time_1 == -999:
                time_1 = tokens[0]
            if time_0 != -999 and time_1 != -999:
                dt = time_1 - time_0
        pgv1 = max(pgv1, abs(tokens[1]))
        pgv2 = max(pgv2, abs(tokens[2]))
        pgv3 = max(pgv3, abs(tokens[3]))
    input_bbp.close()
    
    input_bbp = open(acc_file)
    for line in input_bbp:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        pga1 = max(pga1, abs(tokens[1]))
        pga2 = max(pga2, abs(tokens[2]))
        pga3 = max(pga3, abs(tokens[3]))
    input_bbp.close()
    # Convert cm/s/s to g
    pga1 = pga1 * 0.00101971621
    pga2 = pga2 * 0.00101971621
    pga3 = pga3 * 0.00101971621
    # Calculate frequencies
    nyquist = 1.0 / (2.0 * dt)
    luf = 1.25 * site.high_freq_corner
    huf = 0.8 * min(nyquist, site.low_freq_corner)
    ufb = abs(luf-huf)
    # And periods...
    lup = 1.0 / site.high_freq_corner
    hup = 1.0 / site.low_freq_corner
    upb = abs(lup-hup)
    station["dt"] = dt
    station["num_samples"] = num_samples
    station["nyquist"] = nyquist
    station["luf"] = luf
    station["huf"] = huf
    station["ufb"] = ufb
    station["lup"] = lup
    station["hup"] = hup
    station["upb"] = upb
    station["pga_h1"] = pga1
    station["pga_h2"] = pga2
    station["pga_v"] = pga3
    station["pgv_h1"] = pgv1
    station["pgv_h2"] = pgv2
    station["pgv_v"] = pgv3
    station["rotdnn_fractile"] = 50
    station["damping"] = 5
    station["arias_dur_5_75"] = "TBC"
    station["arias_dur_5_95"] = "TBC"
    station["arias_total"] = "TBC"
    
def collect_rd50_values(station, args):
    """
    Collect RotD50 values for all periods
    """
    rd50_file = os.path.join(args.top_level_outdir, station["rd50_file_name"])
    # Start with an empty list
    rd50_periods = []
    rd50_psa = []
    # Read file
    input_file = open(rd50_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        rd50_periods.append(tokens[0])
        rd50_psa.append(tokens[3])
    # Close file
    input_file.close()

    # All done!
    if "rd50_periods" not in args:
        args.rd50_periods = rd50_periods
    station["rd50"] = rd50_psa

def collect_station_params(site, station, src_files,
                           args, realization):
    """
    Collects parameters for one station
    """
    station["station_name"] = site.scode
    station["recording_station_name"] = "NA"
    station["station_latitude"] = site.lat
    station["station_longitude"] = site.lon
    station["elevation"] = -999
    if site.vs30 is None:
        station["vs30"] = "NA"
        site_class = "NA"
    else:
        station["vs30"] = site.vs30
        if site.vs30 > 1500:
            site_class = "A"
        elif site.vs30 > 760:
            site_class = "B"
        elif site.vs30 > 360:
            site_class = "C"
        elif site.vs30 > 180:
            site_class = "D"
        else:
            site_class = "E"
    station["site_class"] = site_class
    (station["rrup"],
     station["rjb"],
     station["rx"]) = calculate_distances(src_files, site)
    if args.general_method in ["exsim"]:
        station["components"] = 1
    else:
        station["components"] = 3
    station["vel_file_name"] = os.path.join(realization,
                                            "%s.%s.vel.bbp" %
                                            (realization, site.scode))
    station["acc_file_name"] = os.path.join(realization,
                                            "%s.%s.acc.bbp" %
                                            (realization, site.scode))
    station["rd50_file_name"] = os.path.join(realization,
                                             "%s.%s.rd50" %
                                             (realization, site.scode))
    station["h1_azimuth"] = 0
    station["h2_azimuth"] = 90
    station["v_orientation"] = "UP"
    calculate_timeseries_param(station, site, args, realization)
    # Copy files, as needed
    if args.copy_timeseries:
        shutil.copy2(os.path.join(args.top_level_outdir, station["acc_file_name"]),
                     os.path.join(args.output_dir, station["acc_file_name"]))

def collect_realization_params(args, realization):
    """
    Collects parameters for one realization
    """
    indir = os.path.join(args.top_level_indir, realization)
    outdir = os.path.join(args.top_level_outdir, realization)
    src_files = glob.glob("%s/*.src" % (indir))
    stl_file = glob.glob("%s/*.stl" % (indir))[0]
    data = {}

    # Compile data from SRC file(s)
    data["num_src"] = len(src_files)
    # Save info in args too for first realization
    if "num_src" not in args:
        args.num_src = len(src_files)
    for i, src_file in zip(range(1, len(src_files) + 1), src_files):
        src_index = "bbp_src_%d" % (i)
        src_keys = bband_utils.parse_src_file(src_file)
        src_keys["epicenter_latitude"] = "TBC"
        src_keys["epicenter_longitude"] = "TBC"
        src_keys["hypocenter_depth"] = "TBC"
        # Add the GP keys if they are missing
        if "moment_fraction" not in src_keys:
            src_keys["moment_fraction"] = "NA"
        if "max_fault_length" not in src_keys:
            src_keys["max_fault_length"] = "NA"
        if "rupture_delay" not in src_keys:
            src_keys["rupture_delay"] = "NA"
        src_keys["mechanism"] = calculate_mechanism(src_keys["rake"])
        data[src_index] = src_keys

    # Get velocity model data
    html_file = glob.glob("%s/*.html" % (outdir))[0]
    data["vmodel_name"] = get_vmodel_from_html(html_file)
    vel_obj = velocity_models.get_velocity_model_by_name(data["vmodel_name"])
    vmodel_params = vel_obj.get_codebase_params('gp')
    if args.general_method in ["gp", "sdsu", "song"]:
        data["gf_name"] = vmodel_params['GF_NAME']
        data["vs_30"] = "TBC"
    else:
        data["gf_name"] = "NA"
        data["vs_30"] = "NA"

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
        collect_station_params(site, stations[site.scode], src_files,
                               args, realization)
        collect_rd50_values(stations[site.scode], args)
    
    # Save data
    data["stations"] = stations
    
    # Save realization data
    args.data[realization] = data
    
def write_output_data(args):
    """
    This function writes all output to the flat file
    """
    # Output filename
    output_filename = os.path.join(args.output_dir, "bbp-summary-file.csv")
    # Create header
    header = ("bbp_software_version, bbp_simulation_workflow, bbp_site_effects, "
              "record_sequence_number, eq_id, fault_information_source, fault_id, "
              "fault_name, eq_name, eq_location, eq_magnitude, "
              "realization, number_of_src_files")
    for i in range(1, (args.num_src + 1)):
        header = ("%s, src_%d_fault_length, src_%d_fault_width, src_%d_depth_to_top, "
                  "src_%d_strike, src_%d_rake, src_%d_dip, src_%d_latitude, "
                  "src_%d_longitude, src_%d_hypo_along_strike, src_%d_hypo_down_dip, "
                  "src_%d_seed, src_%d_gp_total_fault_length, src_%d_gp_moment_fraction, "
                  "src_%d_gp_timing_info, src_%d_epicenter_latitude, "
                  "src_%d_epicenter_longitude, src_%d_hypocenter_depth, "
                  "src_%d_mechanism" % (header, i, i, i, i,
                                        i, i, i, i, i, i, i,
                                        i, i, i, i, i, i, i))
    header = ("%s, vmodel_name, gf_name, vs30" % (header))
    header = ("%s, station_name, recording_station_name, station_latitude, "
              "station_longitude, station_elevation, station_vs30, "
              "station_class, station_rrup, station_rjb, station_rx, "
              "num_components, acc_file_name, h1_azimuth, h2_azimuth, v_orientation, "
              "dt, num_samples, nyquist, luf, huf, ufb, lup, hup, upb, "
              "pga_h1, pgv_h1, pga_h2, pgv_h2, pga_v, pgv_v, rotdnn_fractile, "
              "damping" % (header))
    for period in args.rd50_periods:
        header = ("%s, %.2f" % (header, period))
    header = ("%s, arias_dur_5_75, arias_dur_5_95, arias_total" % (header))

    # Create first (common) part of the output
    sim_params = ("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s" %
                  (args.bbp_software_info_version,
                   "/".join(args.bbp_software_info_modules),
                   args.bbp_software_info_site,
                   args.general_record_seq_no,
                   args.general_eqid,
                   args.general_fault_information_source,
                   args.general_fault_id,
                   args.general_fault_name,
                   args.general_eq_name,
                   args.general_earthquake_location,
                   str(args.general_magnitude)))

    # Output data
    output_file = open(output_filename, 'w')
    # Print header
    output_file.write("#%s\n" % (header))
    for realization in args.realizations:
        realization_data = args.data[realization]
        station_names = realization_data["station_names"]
        realization_params = ("%s, %d" % (realization, realization_data["num_src"]))
        for i in range(1, (realization_data["num_src"] + 1)):
            src_index = "bbp_src_%d" % (i)
            src_params = realization_data[src_index]
            realization_params = ("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, "
                                  "%s, %s, %s, %s, %s, %s, %s, %s, %s" %
                                  (realization_params, src_params["fault_length"],
                                   src_params["fault_width"],
                                   src_params["depth_to_top"],
                                   src_params["strike"],
                                   src_params["rake"],
                                   src_params["dip"],
                                   src_params["lat_top_center"],
                                   src_params["lon_top_center"],
                                   src_params["hypo_along_stk"],
                                   src_params["hypo_down_dip"],
                                   src_params["seed"],
                                   src_params["max_fault_length"],
                                   src_params["moment_fraction"],
                                   src_params["rupture_delay"],
                                   src_params["epicenter_latitude"],
                                   src_params["epicenter_longitude"],
                                   src_params["hypocenter_depth"],
                                   src_params["mechanism"]))
            realization_params = ("%s, %s, %s, %s" %
                                  (realization_params,
                                   realization_data["vmodel_name"],
                                   realization_data["gf_name"],
                                   realization_data["vs_30"]))
        for station in station_names:
            st_data = realization_data["stations"][station]
            station_params = ("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, "
                              "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, "
                              "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, "
                              "%s, %s" %
                              (station, st_data["recording_station_name"],
                               st_data["station_latitude"],
                               st_data["station_longitude"],
                               st_data["elevation"],
                               st_data["vs30"],
                               st_data["site_class"],
                               st_data["rrup"],
                               st_data["rjb"],
                               st_data["rx"],
                               st_data["components"],
                               st_data["vel_file_name"],
                               st_data["h1_azimuth"],
                               st_data["h2_azimuth"],
                               st_data["v_orientation"],
                               st_data["dt"],
                               st_data["num_samples"],
                               st_data["nyquist"],
                               st_data["luf"],
                               st_data["huf"],
                               st_data["ufb"],
                               st_data["lup"],
                               st_data["hup"],
                               st_data["upb"],
                               st_data["pga_h1"],
                               st_data["pgv_h1"],
                               st_data["pga_h2"],
                               st_data["pgv_h2"],
                               st_data["pga_v"],
                               st_data["pgv_v"],
                               st_data["rotdnn_fractile"],
                               st_data["damping"]))
            for period in st_data["rd50"]:
                station_params = ("%s, %s" % (station_params, period))
            station_params = ("%s, %s, %s, %s" % (station_params,
                                                  st_data["arias_dur_5_75"],
                                                  st_data["arias_dur_5_95"],
                                                  st_data["arias_total"]))
            output_file.write("%s, %s, %s\n" %
                              (sim_params, realization_params, station_params))
    # All done
    output_file.close()
    
def create_flat_file_from_cluster():
    """
    Create a flat file from a cluster simulation
    """
    # Get all we need from the command-line
    args = parse_arguments()

    # Figure out top-level directories
    args.top_level_indir = os.path.join(args.input_dir, "Sims", "indata")
    args.top_level_outdir = os.path.join(args.input_dir, "Sims", "outdata")
    args.realizations = sorted(os.listdir(args.top_level_indir))
    args.data = {}

    # Create top-level output directory
    bband_utils.mkdirs([args.output_dir], print_cmd=False)
    
    # Collect simulation-wide parameters
    collect_simulation_params(args)

    # Collect parameters for each realization
    for realization in args.realizations:
        print("==> Processing realization: %s..." % (realization))
        # Create output directory for realization if requested to copy seismograms
        if args.copy_timeseries:
            bband_utils.mkdirs([os.path.join(args.output_dir, realization)],
                               print_cmd=False)
        collect_realization_params(args, realization)
    
    # Write flat file
    write_output_data(args)

if __name__ == '__main__':
    create_flat_file_from_cluster()
