#!/bin/env python
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

This module plots a Broadband velocity model profile
"""
from __future__ import division, print_function

# Import Python Modules
import os
import sys
import math
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def create_velocity_plot(data, output_file,
                         ylabel=None, xlabel=None,
                         title=None, legend=None):
    """
    Creates the velocity model plot
    """
    # Plot parameters
    width = 7
    height = 10
    vplabel = "Vp (km/s)"
    vslabel = "Vs (km/s)"
    densitylabel = "Density (g/cm^3)"
    vpcolor = "r"
    vscolor = "b"
    densitycolor = "g"
    startingdepth = 0.0
    todepth = 60000.0

    # Data
    thicklist = data[0]
    vplist = data[1]
    vslist = data[2]
    rholist = data[3]

    figure = plt.figure(figsize=(width, height), dpi=100)
    if ylabel != None:
        plt.ylabel(ylabel, fontsize=14)

    if xlabel != None:
        plt.xlabel(xlabel, fontsize=14)

    if title != None:
        plt.title(title)

    if legend != None:
        plt.legend(legend, loc='lower left')

    # Calculate max_x value
    max_x = 0
    max_x = max(max_x, max(vplist))
    max_x = max(max_x, max(vslist))
    max_x = max(max_x, max(rholist))
    max_x = math.ceil(max_x)

    # Go over data to create points to plot
    new_thick_list = [0.0]
    new_vp_list = [0.0]
    new_vs_list = [0.0]
    new_rho_list = [0.0]
    for thick, vp, vs, rho in zip(thicklist, vplist, vslist, rholist):
        # Convert thick from km to m
        thick = thick * 1000
        # Figure out current thickness
        curr_thick = new_thick_list[-1] + thick
        new_thick_list.append(new_thick_list[-1])
        new_vp_list.append(vp)
        new_vs_list.append(vs)
        new_rho_list.append(rho)
        new_thick_list.append(curr_thick)
        new_vp_list.append(vp)
        new_vs_list.append(vs)
        new_rho_list.append(rho)

    # Create plots
    figure.add_subplot(1, 1, 1).plot(new_vp_list, new_thick_list, '-',
                                     color=vpcolor, label=vplabel)
    figure.add_subplot(1, 1, 1).plot(new_vs_list, new_thick_list, '-',
                                     color=vscolor, label=vslabel)
    figure.add_subplot(1, 1, 1).plot(new_rho_list, new_thick_list, '-',
                                     color=densitycolor, label=densitylabel)
    # Add legend box
    plt.legend(loc="lower left")

    # Set axis limits
    if plt.ylim()[0] < plt.ylim()[1]:
        plt.gca().invert_yaxis()

    if max_x > plt.xlim()[1]:
        plt.xlim(0, math.ceil(max_x / 0.5) * 0.5)

    plt.axis([0, max_x, int(todepth), int(startingdepth)])

    # All done, save plot
    plt.savefig(output_file)

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Creates a map plot with "
                                     " a fault and stations.")
    parser.add_argument("-o", "--output", dest="outfile", required=True,
                        help="output png file")
    parser.add_argument("-v", "--velocity-model",
                        dest="velocity_model", required=True,
                        help="velocity model file to plot")
    parser.add_argument("--title", "-t",
                        dest="plot_title",
                        help="title for plot")
    args = parser.parse_args()

    return args

def read_data(input_file):
    """
    Reads the velocity model specified in input_file
    """
    thick = []
    vp = []
    vs = []
    rho =[]
    vel_model = open(input_file, 'r')
    for line in vel_model:
        line = line.strip()
        # Skip blank lines
        if not line:
            continue
        tokens = line.split()
        # Only read lines that have all 6 tokens
        if len(tokens) != 6:
            continue
        tokens = [float(token) for token in tokens]
        thick.append(tokens[0])
        vp.append(tokens[1])
        vs.append(tokens[2])
        rho.append(tokens[3])
    # All done
    vel_model.close()
    data = [thick, vp, vs, rho]
    return data

def plot_velocity_main():
    """
    Main function for generating a velocity model plot
    """
    # Parse command-line options
    args = parse_arguments()
    data = read_data(args.velocity_model)

    create_velocity_plot(data=data, output_file=args.outfile,
                         ylabel="Depth (m)",
                         xlabel="Units (see legend)",
                         title=args.plot_title,
                         legend=None)

# ============================ MAIN ==============================
if __name__ == "__main__":
    plot_velocity_main()
# end of main program
