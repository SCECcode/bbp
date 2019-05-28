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

Program used to learn how many realizations are needed
for a method to converge.

This program is based on the original bbp_converge script
from Karen Assatourians, Western University
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import math
import curses
import random
import shutil
import argparse
import numpy as np
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import pylab

# Import Broadband modules
import bband_utils
import plot_config

DEFAULT_SEED = 123456789

def update_progress_bar(progress, total):
    """
    Keeps the progress bar moving
    """
    bar_size = 40
    completed = int(progress * bar_size / total)
    missing = bar_size - completed
    progress_bar = "%s%s" % ("*" * completed, "-" * missing)
    print("\r Iteration: %s : %d" % (progress_bar, progress), end="")
    sys.stdout.flush()

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Learn how many "
                                     " realizations are needed "
                                     " for methods to converge.")
    parser.add_argument("--input_dir", "-i", dest="input_dir",
                        required=True, help="input directory")
    parser.add_argument("-o", "--output",
                        dest="output_file", required=True,
                        help="output png file")
    parser.add_argument("--limit", "-l", type=float,
                        dest="limit", default=0.02,
                        help="difference limit")
    parser.add_argument("--ns", type=int, default=10000,
                        dest="sampling", help="number of sampling")
    parser.add_argument("-c", "--codebase", required=True,
                        dest="codebase",
                        help="method used for the simulation")
    parser.add_argument("--colormap", default="Paired",
                        dest="colormap",
                        help="matplotlib colormap to use")
    args = parser.parse_args()

    return args

def read_input_bias_data(input_dir):
    """
    Read the bias data from all realizations
    """
    periods = []
    data = []
    event_label = None

    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization)
        bias_file = glob.glob("%s%s*-rotd50.bias" % (basedir, os.sep))
        if len(bias_file) != 1:
            raise bband_utils.ProcessingError("Bias file not found for "
                                              "realization %s!" % (realization))
        bias_file = bias_file[0]
        # Let's capture the event label
        if event_label is None:
            event_label = os.path.basename(bias_file).split("-")[0]

        input_file = open(bias_file, 'r')
        cur_periods = []
        cur_data = []
        for line in input_file:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith("#") or line.startswith("%") or not line:
                continue
            tokens = [float(token) for token in line.split()]
            cur_periods.append(tokens[0])
            cur_data.append(tokens[1])
        # Close input_file
        input_file.close()
        # Keep list of periods if not already done
        if not periods:
            periods = cur_periods
        # And keep data
        data.append(cur_data)

    bias_data = {}
    bias_data["num_periods"] = len(periods)
    bias_data["periods"] = periods
    bias_data["num_realizations"] = len(realizations)
    bias_data["data"] = data
    bias_data["event_label"] = event_label

    return bias_data

def find_gavg(bias_data):
    """
    Calculate averages
    """
    gavg = np.zeros(bias_data["num_periods"])
    data = bias_data["data"]
    num_realizations = float(bias_data["num_realizations"])
    for realization in data:
        for cur_period in range(0, bias_data["num_periods"]):
            gavg[cur_period] = (gavg[cur_period] +
                                realization[cur_period] / num_realizations)

    bias_data["gavg"] = gavg
    return bias_data

def calculate_probabilities(bias_data, limit, sampling):
    """
    Calculate the probabilities
    """
    ratios_by_period = [[] for _ in range(0, bias_data["num_periods"])]

    random.seed(DEFAULT_SEED)
    for cur_realization in range(1, bias_data["num_realizations"] + 1):
        update_progress_bar(cur_realization, bias_data["num_realizations"])
        ratio = np.zeros(bias_data["num_periods"])
        for _ in range(0, sampling):
            avg = np.zeros(bias_data["num_periods"])
            for _ in range(0, cur_realization):
                random_int = (int(random.random() * bias_data["num_realizations"]))
                if random_int > bias_data["num_realizations"]:
                    random_int = bias_data["num_realizations"]
                for nper in range(0, bias_data["num_periods"]):
                    avg[nper] = (avg[nper] +
                                 bias_data["data"][random_int][nper] / float(cur_realization))
            for nper in range(0, bias_data["num_periods"]):
                if abs(avg[nper] - bias_data["gavg"][nper]) <= limit:
                    ratio[nper] = ratio[nper] + 1
        ratio = ratio / sampling

        # Now re-shuffle so that we split by periods not realizations
        for val, dest in zip(ratio, ratios_by_period):
            dest.append(val)

    # Save calculated data
    bias_data["ratios"] = ratios_by_period
    print()

    return bias_data

def plot_results(bias_data, codebase, colormap, output_file):
    """
    Generate plot showing results calculated from all realizations
    """
    fig, ax = pylab.plt.subplots()
    xlocs = range(1, bias_data["num_periods"] + 1)
    ylocs = list(np.full(bias_data["num_periods"],
                         bias_data["num_realizations"]))
    bars = ax.bar(xlocs, ylocs)
    ax = bars[0].axes
    lim = ax.get_xlim() + ax.get_ylim()
    for bar, ratio in zip(bars, bias_data["ratios"]):
        rev_ratio = ratio[::-1]
        gradient = np.atleast_2d(rev_ratio).T
        bar.set_zorder(1)
        bar.set_facecolor("none")
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        im = ax.imshow(gradient, extent=[x, x + w, y, y + h],
                       cmap=pylab.get_cmap(colormap),
                       vmin=0.0, vmax=1.0,
                       aspect="auto", zorder=0)
    ax.axis(lim)
    fig.colorbar(im)
    pylab.xlim(0, bias_data["num_periods"] + 2)
    pylab.plt.xticks([val + 0.5 for val in xlocs],
                     bias_data["periods"],
                     rotation='vertical', fontsize=6)
    #frame1 = pylab.gca()
    #frame1.axes.xaxis.set_ticklabels([])
    #frame1.axes.yaxis.set_ticklabels([])
    pylab.xlabel("Periods")
    pylab.ylabel("Number of Realizations")
    pylab.title(("Convergence Plot - %s Method - %s" %
                 (codebase, bias_data["event_label"])),
                size=10)

    # Save plot
    fig.savefig(output_file, format="png",
                transparent=False, dpi=plot_config.dpi)

def main():
    """
    Figure out how many realizations are needed for methods to
    converge
    """
    # Parse command-line options
    args = parse_arguments()

    # Copy options
    limit = args.limit
    sampling = args.sampling
    base_input_dir = args.input_dir
    output_file = args.output_file
    colormap = args.colormap
    codebase = args.codebase

    # Read input data
    input_dir = os.path.join(base_input_dir, "Sims", "outdata")
    bias_data = read_input_bias_data(input_dir)
    bias_data = find_gavg(bias_data)
    bias_data = calculate_probabilities(bias_data, limit, sampling)
    plot_results(bias_data, codebase, colormap, output_file)

if __name__ == "__main__":
    main()
