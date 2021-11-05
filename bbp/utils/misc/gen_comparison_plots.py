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

Program to plot seismogram comparison between two BBP files
"""
# Import Python modules
from __future__ import division, print_function
import os
import sys
import glob
import argparse

import matplotlib as mpl
mpl.use("AGG", warn=False)
import pylab
import numpy as np
from scipy.integrate import cumtrapz

# Import BBP modules
import plot_config
from station_list import StationList

def integrate(data, dt):
    """
    Compute integral of a numpy array
    initial condition assumed 0
    result has same size as input
    """
    data = np.array(data)
    newdata = list(cumtrapz(data, dx=dt, initial=0) + data[0] * dt / 2.0)
    return newdata

def derivative(data, dt):
    """
    Compute derivative of a timeseries
    """
    data = np.array(data)
    newdata = np.insert(data, 0, 0)
    newdata = list(np.diff(newdata) / dt)
    return newdata

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Generates timeseries plots "
                                     "of velocity BBP files.")
    parser.add_argument("--list", "--station-list", dest="station_list",
                        help="station list containing list of files to process")
    parser.add_argument("--xmin", dest="xmin",
                        help="minimum x value (default: 0)")
    parser.add_argument("--xmax", dest="xmax",
                        help="maximum x value (default: entire timeseries)")
    parser.add_argument("--output-dir", "-o", dest="output_dir",
                        help="output directory for the plots")
    parser.add_argument("--dpi", dest="dpi",
                        help="dpi to use in the plots "
                        "(default: dpi in plot_config)")
    parser.add_argument("-c", "--method", dest="method",
                        help="specify the method used in the simulations")
    parser.add_argument("--label1", dest="label1",
                        help="label for first dataset")
    parser.add_argument("--label2", dest="label2",
                        help="label for second dataset")
    parser.add_argument("--input", dest="input",
                        help="input data files: dis, vel, acc")
    parser.add_argument("--dis", action="store_true", dest="dis",
                        help="generate displacement plots")
    parser.add_argument("--vel", action="store_true", dest="vel",
                        help="generate velocity plots")
    parser.add_argument("--acc", action="store_true", dest="acc",
                        help="generate acceleration plots")
    parser.add_argument('input_dirs', nargs='*')
    args = parser.parse_args()

    # dirs = args.input_dirs
    if len(args.input_dirs) < 1 or len(args.input_dirs) > 2:
        print("ERROR: Please provide one or two directories for input files")
        sys.exit(-1)

    if args.station_list is None:
        print("ERROR: Please provide station list file")
        sys.exit(-1)

    if args.input is None:
        print("ERROR: Please provide input type")
        sys.exit(-1)

    if not (args.dis or args.vel or args.acc):
        print("ERROR: Please select plot type to generate")
        sys.exit(-1)

    if (args.dis and args.vel) or (args.dis and args.acc) or (args.vel and args.acc):
        print("ERROR: Please select only one plot type to generate")
        sys.exit(-1)

    params = {}
    params["input_dirs"] = args.input_dirs
    params["station_list"] = args.station_list
    if args.xmin is None:
        params["xmin"] = 0
    else:
        params["xmin"] = float(args.xmin)
    if args.xmax is None:
        params["xmax"] = None
    else:
        params["xmax"] = float(args.xmax)
    if args.output_dir is None:
        params["output_dir"] = "."
    else:
        params["output_dir"] = args.output_dir
    if args.dpi is None:
        params["dpi"] = plot_config.dpi
    else:
        params["dpi"] = args.dpi
    if args.method is None:
        params["method"] = None
    else:
        params["method"] = args.method

    if args.input.lower() == "dis":
        params["extension"] = ".dis.bbp"
        params["input"] = "dis"
    elif args.input.lower() == "vel":
        params["extension"] = ".vel.bbp"
        params["input"] = "vel"
    elif args.input.lower() == "acc":
        params["extension"] = ".acc.bbp"
        params["input"] = "acc"
    else:
        print("ERROR: Invalid input data type, please specify dis, vel or acc")
        sys.exit(-1)

    if args.dis == True:
        params["output"] = "dis"
        params["units"] = "(cm)"
    if args.vel == True:
        params["output"] = "vel"
        params["units"] = "(cm/s)"
    if args.acc == True:
        params["output"] = "acc"
        params["units"] = "(cm/s/s)"

    if args.label1 is None:
        params["label1"] = "Set 1"
    else:
        params["label1"] = args.label1
    if args.label2 is None:
        params["label2"] = "Set 2"
    else:
        params["label2"] = args.label2

    return params

def read_seismogram_file(filename, params):
    """
    This function reads a seismogram and returns 4 lists with the
    horizontal components (ns and ew), vertical, and the timestamps
    """
    # Start empty
    time_comp = []
    ns_comp = []
    ew_comp = []
    ver_comp = []

    # Read file
    seis_file = open(filename, 'r')
    for line in seis_file:
        # Remove leading spaces
        line = line.strip()
        # Skip comments
        if line.startswith('#') or line.startswith('%'):
            continue
        tmp = line.split()
        if len(tmp) < 4:
            print("Error reading seismogram in file %s" % (filename))
            sys.exit(1)
        time_comp.append(float(tmp[0]))
        ns_comp.append(float(tmp[1]))
        ew_comp.append(float(tmp[2]))
        ver_comp.append(float(tmp[3]))
    # Close file
    seis_file.close()

    if params["input"] == params["output"]:
        # All done
        return (time_comp, ns_comp, ew_comp, ver_comp)

    # Change units as needed
    delta_t = time_comp[1] - time_comp[0]

    # Set input variables
    if params["input"] == "dis":
        dis_ns = ns_comp
        dis_ew = ew_comp
        dis_ud = ver_comp
    elif params["input"] == "vel":
        vel_ns = ns_comp
        vel_ew = ew_comp
        vel_ud = ver_comp
    elif params["input"] == "acc":
        acc_ns = ns_comp
        acc_ew = ew_comp
        acc_ud = ver_comp
    else:
        print("ERROR: input type must be one of: dis, vel, acc")
        sys.exit(-1)

    if params["input"] == "dis":
        vel_ns = derivative(dis_ns, delta_t)
        vel_ew = derivative(dis_ew, delta_t)
        vel_ud = derivative(dis_ud, delta_t)
    if params["output"] == "acc":
        acc_ns = derivative(vel_ns, delta_t)
        acc_ew = derivative(vel_ew, delta_t)
        acc_ud = derivative(vel_ud, delta_t)
    if params["input"] == "acc":
        vel_ns = integrate(acc_ns, delta_t)
        vel_ew = integrate(acc_ew, delta_t)
        vel_ud = integrate(acc_ud, delta_t)
    if params["output"] == "dis":
        dis_ns = integrate(vel_ns, delta_t)
        dis_ew = integrate(vel_ew, delta_t)
        dis_ud = integrate(vel_ud, delta_t)

    if params["output"] == "dis":
        return (time_comp, dis_ns, dis_ew, dis_ud)
    elif params["output"] == "vel":
        return (time_comp, vel_ns, vel_ew, vel_ud)
    elif params["output"] == "acc":
        return (time_comp, acc_ns, acc_ew, acc_ud)

    print("ERROR: output type must be one of: dis, vel, acc")
    sys.exit(-1)

def plot_comparison(files, outfile, params):
    """
    Creates a comparison plot between a set of files.
    Outputs the plot to outfile
    """
    blue = '#3182bd'
    red = '#ef3b2c'
    alpha = 0.7
    alpha_ = 0.85

    ts1, ns1, ew1, ud1 = read_seismogram_file(files[0], params)
    # Read second set, if provided
    if len(files) > 1:
        ts2, ns2, ew2, ud2 = read_seismogram_file(files[1], params)

    if params["xmax"] is None:
        if len(files) > 1:
            xmax = max(max(ts1), max(ts2))
        else:
            xmax = max(ts1)
    else:
        xmax = params["xmax"]

    _, (ax1, ax2, ax3) = pylab.subplots(3, sharex=True)
    if len(files) > 1:
        if np.max(np.abs(ns1)) >= np.max(np.abs(ns2)):
            ax1.plot(ts1, ns1, c=blue, label=params["label1"])
            ax1.plot(ts1, ns2, c=red, label=params["label2"], alpha=alpha)
        else:
            ax1.plot(ts1, ns2, c=red, label=params["label2"], alpha=alpha_)
            ax1.plot(ts1, ns1, c=blue, label=params["label1"], alpha=alpha_)
        if np.max(np.abs(ew1)) >= np.max(np.abs(ew2)):
            ax2.plot(ts1, ew1, c=blue, label=params["label1"])
            ax2.plot(ts1, ew2, c=red, label=params["label2"], alpha=alpha)
        else:
            ax2.plot(ts1, ew2, c=red, label=params["label2"], alpha=alpha_)
            ax2.plot(ts1, ew1, c=blue, label=params["label1"], alpha=alpha_)
        if np.max(np.abs(ud1)) >= np.max(np.abs(ud2)):
            ax3.plot(ts1, ud1, c=blue, label=params["label1"])
            ax3.plot(ts1, ud2, c=red, label=params["label2"], alpha=alpha)
        else:
            ax3.plot(ts1, ud2, c=red, label=params["label2"], alpha=alpha_)
            ax3.plot(ts1, ud1, c=blue, label=params["label1"], alpha=alpha_)
    else:
        ax1.plot(ts1, ns1, c=blue, label=params["label1"])
        ax2.plot(ts1, ew1, c=blue, label=params["label1"])
        ax3.plot(ts1, ud1, c=blue, label=params["label1"])
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    ax3.legend(loc='upper right')
    ax1.grid(ls=":")
    ax2.grid(ls=":")
    ax3.grid(ls=":")
    ax1.set_ylabel("N/S %s" % (params["units"]))
    ax2.set_ylabel("E/W %s" % (params["units"]))
    ax3.set_ylabel("U/D %s" % (params["units"]))
    if params["method"] is None:
        ax1.set_title("Station: %s - Vs30: %d" % (params["station"],
                                                  int(params["vs30"])))
    else:
        ax1.set_title("Station: %s - Vs30: %d - Method: %s" %
                      (params["station"], int(params["vs30"]), params["method"]))
    ax3.set_xlabel("Time (s)")

    # Set xmin, xmax
    for cur_x in (ax1, ax2, ax3):
        cur_x.set_xlim(params["xmin"], xmax)

    # Save plot
    pylab.gcf().set_size_inches(18, 12)
    pylab.savefig(outfile, format="png", dpi=params["dpi"])

    # Remember to close plot
    pylab.close()

def main():
    """
    Main function, parses command-line options, calls plotting code
    """
    params = parse_arguments()

    slo = StationList(params["station_list"])
    station_list = slo.getStationList()

    for site in station_list:
        station = site.scode
        print("==> Processing data for station: %s" % (station))

        filename = "%s%s" % (station, params["extension"])
        files = []
        for input_dir in params["input_dirs"]:
            input_file = glob.glob("%s%s*%s" % (input_dir, os.sep, filename))
            if len(input_file) == 1:
                input_file = input_file[0]
            else:
                print("[Error]: Cannot find data file for station %s" % (station))
                sys.exit(-1)

            files.append(input_file)

        outfile = os.path.join(params["output_dir"],
                               "%s-%s-comparison.png" %
                               (station, params["output"]))
        params["station"] = station
        params["vs30"] = site.vs30

        plot_comparison(files, outfile, params)

if __name__ == "__main__":
    main()
