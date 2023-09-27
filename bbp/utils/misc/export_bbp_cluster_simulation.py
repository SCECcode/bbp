#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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
import datetime
import numpy as np
from scipy import integrate

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
CM2G = 980.664999

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
    parser.add_argument("--prefix", "-p", dest="prefix",
                        default="BBP_Study",
                        help="prefix for output files")
    parser.add_argument("--suffix", "-s", dest="suffix",
                        default="",
                        help="suffix for output files")
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

    args.general_eqid = "-999"

def calculate_vs30(vmodel_file):
    """
    Calculates the Vs30 from the velocity model file
    """
    # Need to calculate Vs30 by adding each layer's up to 30m
    total_time = 0.0
    remaining_width = 30.0

    input_file = open(vmodel_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        tokens = [float(item) for item in line.split()]
        if len(tokens) != 6:
            continue
        layer_width = tokens[0] * 1000 # Convert to meters
        layer_vs = tokens[2] * 1000
        if layer_width <= remaining_width:
            remaining_width = remaining_width - layer_width
            total_time = total_time + (layer_width / layer_vs)
        else:
            total_time = total_time + (remaining_width / layer_vs)
            remaining_width = 0.0
        # Check if all done
        if remaining_width == 0.0:
            break
    input_file.close()

    # Calculate Vs30 based on total_time
    return (30.0 / total_time)

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

def read_bbp(bbp_file):
    """
    This function reads the input bbp_file and returns 4 arrays
    containing the timestamps and 3 time series. This function
    converts all BBP files from cm/2^2 to g
    """
    times = []
    comp1 = []
    comp2 = []
    comp3 = []

    ifile = open(bbp_file, 'r')
    for line in ifile:
        line = line.strip()
        # Skip comments
        if line.startswith('%') or line.startswith('#'):
            continue
        pieces = [float(x) for x in line.split()]
        times.append(pieces[0])
        comp1.append(pieces[1] / CM2G)
        comp2.append(pieces[2] / CM2G)
        comp3.append(pieces[3] / CM2G)

    # Close input file
    ifile.close()

    # All done, return arrays
    return times, comp1, comp2, comp3

def calculate_arias(F, dt, percent):
    """
    For a given motion, this function will tell you at what time a
    given percentage of arias intensity is reached (if time starts
    at 0 sec)
    """
    n = len(F)

    a_i = [pow(value, 2) for value in F]
    I = integrate.cumtrapz(a_i) * dt
    # Arias Intensity
    Ia = (F[0]**2) * dt / 2.0 + I[n-2] + (F[n-1]**2) * dt / 2.0
    It = (percent / 100.0) * Ia

    if I[0] < It:
        index = len(I) - len(I[I >= It])
        if index == len(I):
            index = index - 1
    else:
        index = 0

    t = index * dt
    return t, index, It

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
    # Calculate arias duration values
    bbp_data = read_bbp(acc_file)
    t5_h1, _, _ = calculate_arias(bbp_data[1], dt, 5)
    t5_h2, _, _ = calculate_arias(bbp_data[2], dt, 5)
    t5_v, _, _ = calculate_arias(bbp_data[3], dt, 5)
    t75_h1, _, _ = calculate_arias(bbp_data[1], dt, 75)
    t75_h2, _, _ = calculate_arias(bbp_data[2], dt, 75)
    t75_v, _, _ = calculate_arias(bbp_data[3], dt, 75)
    t95_h1, _, _ = calculate_arias(bbp_data[1], dt, 95)
    t95_h2, _, _ = calculate_arias(bbp_data[2], dt, 95)
    t95_v, _, _ = calculate_arias(bbp_data[3], dt, 95)
    # Calculate times
    station["ai_h1"] = -999
    station["ai_h2"] = -999
    station["ai_v"] = -999
    station["ad5_75_h1"] = t75_h1 - t5_h1
    station["ad5_75_h2"] = t75_h2 - t5_h2
    station["ad5_75_v"] = t75_v - t5_v
    station["ad5_95_h1"] = t95_h1 - t5_h1
    station["ad5_95_h2"] = t95_h2 - t5_h2
    station["ad5_95_v"] = t95_v - t5_v

    station["time_series_dt"] = dt
    station["time_series_num_samples"] = num_samples
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
    station["rotdnn_fractile"] = "PSA_RotD50"
    station["damping"] = 0.05
    station["arias_dur_5_75"] = "-999"
    station["arias_dur_5_95"] = "-999"
    station["arias_total"] = "-999"

def collect_rd50_values(station, args):
    """
    Collect RotD50 values for all periods
    """
    rd50_file = os.path.join(args.top_level_outdir, station["rd50_file_name"])
    rd50_vertical_file = os.path.join(args.top_level_outdir,
                                      station["rd50_vertical_file_name"])

    # Start with an empty list
    rd50_periods = []
    rd50_psa_h1 = []
    rd50_psa_h2 = []
    rd50_psa_v = []
    rd50_psa_rd50 = []

    # Read horizontal psa file
    input_file = open(rd50_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        rd50_periods.append(tokens[0])
        rd50_psa_h1.append(tokens[1])
        rd50_psa_h2.append(tokens[2])
        rd50_psa_rd50.append(tokens[3])
    # Close file
    input_file.close()

    # Read vertical psa file
    input_file = open(rd50_vertical_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        rd50_psa_v.append(tokens[1])
    # Close file
    input_file.close()

    # All done!
    if "rd50_periods" not in args:
        args.rd50_periods = rd50_periods
    station["psa_h1"] = rd50_psa_h1
    station["psa_h2"] = rd50_psa_h2
    station["psa_v"] = rd50_psa_v
    station["rd50"] = rd50_psa_rd50

def collect_rd100_values(station, args):
    """
    Collect RotD100 values for all periods
    """
    rd100_file = os.path.join(args.top_level_outdir, station["rd100_file_name"])

    # Skip if RD100 file doesn't exist
    if not os.path.isfile(rd100_file):
        # RotD100 file not available
        station["rd100"] = None
        return

    # Start with an empty list
    rd100_psa_rd100 = []
    # Read file
    input_file = open(rd100_file, 'r')
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        tokens = [float(token) for token in line.split()]
        rd100_psa_rd100.append(tokens[4])
    # Close file
    input_file.close()

    # All done!
    station["rd100"] = rd100_psa_rd100

def collect_station_params(site, station, src_files,
                           args, realization, vs_30):
    """
    Collects parameters for one station
    """
    station["sim_station_name"] = site.scode
    station["sim_station_latitude"] = site.lat
    station["sim_station_longitude"] = site.lon
    station["sim_station_elevation"] = -999.0
    if isinstance(vs_30, int):
        station["target_station_vs30"] = vs_30
        if vs_30 > 1500:
            site_class = "A"
        elif vs_30 > 760:
            site_class = "B"
        elif vs_30 > 360:
            site_class = "C"
        elif vs_30 > 180:
            site_class = "D"
        else:
            site_class = "E"
    else:
        station["target_station_vs30"] = "-888"
        site_class = "-888"

    station["target_station_nehrp_class"] = site_class
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
    station["rd50_vertical_file_name"] = os.path.join(realization,
                                                      "%s.%s.rd50.vertical" %
                                                      (realization, site.scode))
    station["rd100_file_name"] = os.path.join(realization,
                                              "%s.%s.rd100" %
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
    html_file = glob.glob("%s/*.html" % (outdir))[0]
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
    site_list = slo.get_station_list()
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

def write_output_data(args):
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
    output_psa_period_table_filename = os.path.join(args.output_dir,
                                                    "%s_%s%s_PSA_Period_Table%s" %
                                                    (output_filename_prefix,
                                                     output_filename_date,
                                                     output_filename_suffix,
                                                     output_filename_extension))

    # Create header
    header = ("acc_file_name,bbp_software_version,sim_simulation_workflow,"
              "sim_method_short_name,sim_site_effects,"
              "eq_id,eq_magnitude,"
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
              "num_components,h1_azimuth,h2_azimuth,v_orientation,"
              "dt,num_samples,nyquist,luf,huf,ufb,lup,hup,upb,"
              "pga_h1,pga_h2,pga_v,pgv_h1,pgv_h2,pgv_v" % (header))
    header = ("%s,ai_h1,ai_h2,ai_v,aid5_75_h1,aid5_75_h2,aid5_75_v,"
              "aid5_95_h1,aid5_95_h2,aid5_95_v" % (header))

    header_psa = "acc_file_name,intensity_measure,damping"
    header_periods = "T%dp%03d" % (int(args.rd50_periods[0]),
                                   args.rd50_periods[0] % 1 * 1000)
    for period in args.rd50_periods[1:]:
        header_periods = ("%s,T%dp%03d" % (header_periods,
                                           int(period),
                                           (period % 1 * 1000)))
    header_psa = "%s,%s" % (header_psa, header_periods)

    # Create first (common) part of the output
    sim_params = ('"%s","%s","%s","%s","%s",%s' %
                  (args.bbp_software_info_version,
                   "/".join(args.bbp_software_info_modules),
                   args.general_method,
                   args.bbp_software_info_site,
                   args.general_eqid,
                   str(args.general_magnitude)))

    # Output PSA period table
    output_file = open(output_psa_period_table_filename, 'w')
    output_file.write("%s\n" % (header_periods))
    output_file.write("%.3f" % (args.rd50_periods[0]))
    for period in args.rd50_periods[1:]:
        output_file.write(",%.3f" % (period))
    output_file.write("\n")
    output_file.close()

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
                              '%s,%s,"%s",%s,%s,%s,%s,%s,%s,'
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
            output_file.write('"%s",%s,%s,%s\n' %
                              (st_data["acc_file_name"], sim_params,
                               realization_params, station_params))

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
