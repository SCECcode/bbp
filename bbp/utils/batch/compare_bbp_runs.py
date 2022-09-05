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

Tool to compare bias plots from Broadband simulations done with the
bbp parallel scripts on a cluster
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import operator
import argparse

def compute_avg_bias(bfn):
    """
    Compute the average bias for the bias_file
    """
    bias = 0.0
    points = 0

    bias_file = open(bfn)
    for line in bias_file:
        values = line.split()
        # Each line should have 2 items
        if len(values) != 2:
            continue
        val = float(values[1])
        # Skip added zeroes
        if val == 0.0:
            continue
        bias = bias + abs(val)
        points = points + 1

    # Done adding, calculate the average
    if points > 0:
        bias = bias / points
    # Close file
    bias_file.close()

    # Return bias
    return bias

def print_results_psa(results_psa):
    """
    Print PSA results
    """
    # Print the results
    print("Rank    Simulation    Average Bias(PSA)")
    rank = 0
    for key, val in sorted(results_psa.items(),
                           key=operator.itemgetter(1)):
        rank = rank + 1
        print("%4d    %10s        %8f" % (rank, key, val))

def print_results_fas(results_fas):
    """
    Print PSA and FAS results
    """
    # Print the results
    print("Rank    Simulation    Average Bias(FAS)")
    rank = 0
    for key, val in sorted(results_fas.items(),
                           key=operator.itemgetter(1)):
        rank = rank + 1
        print("%4d    %10s        %8f" % (rank, key, val))

def print_results_combined(results_psa, results_fas):
    """
    Create comparison including PSA and FAS results
    """
    sorted_realizations = sorted(results_fas)
    sorted_psa = sorted(results_psa.items(),
                        key=operator.itemgetter(1))
    sorted_fas = sorted(results_fas.items(),
                        key=operator.itemgetter(1))
    rank_psa = list(dict(sorted_psa).keys())
    rank_fas = list(dict(sorted_fas).keys())

    print("Simulation  Rank(PSA)  Rank(FAS)   Average Bias(PSA)   Average Bias(FAS)")
    for realization in sorted_realizations:
        print("%10s   %6d    %6d     %8f       %8f" % (realization,
                                                       rank_psa.index(realization) + 1,
                                                       rank_fas.index(realization) + 1,
                                                       results_psa[realization],
                                                       results_fas[realization]))
        
def compare_runs(top_dir, args, output_file=None):
    """
    For all sub-directories in top_dir, look for the output bias
    results, calculate stats, and rank it. The result goes to
    output_file, or to stdout if output_file is not provided.
    """
    # Start with empty lists
    results_psa = {}
    results_fas = {}
    fas = args.fas
    psa_fas_comp = args.comp

    for item in os.listdir(top_dir):
        sim_dir = os.path.join(top_dir, item)

        # Only look at directories
        if not os.path.isdir(sim_dir):
            continue

        files_psa = glob.glob(os.path.join(sim_dir, "*rotd50.bias"))
        files_fas = glob.glob(os.path.join(sim_dir, "*seas.bias"))
        if len(files_psa) != 1:
            # PSA bias file not found
            continue
        if fas and len(files_fas) != 1:
            # FAS bias file not found
            continue

        # Calculate stats for this simulation and add to results
        results_psa[item] = compute_avg_bias(files_psa[0])
        if fas:
            results_fas[item] = compute_avg_bias(files_fas[0])

    if not fas:
        print_results_psa(results_psa)
    else:
        if psa_fas_comp:
            # Combine FAS and PSA results
            print_results_combined(results_psa, results_fas)
        else:
            print_results_fas(results_fas)

def parse_arguments():
    """
    This function takes care of parsing the command-line arguments and
    asking the user for any missing parameters that we need
    """
    parser = argparse.ArgumentParser(description="Compare BBP "
                                     " realizations by reading "
                                     " PSA and FAS bias files.")
    parser.add_argument("--dir", "-d", dest="input_dir",
                        required=True, help="input directory")
    parser.add_argument("-o", "--output-dir", required=False,
                        dest="output_dir", default=None,
                        help="output directory")
    parser.add_argument("--fas", action='store_true',
                        help="Add FAS comparison")
    parser.add_argument("--comp", action="store_true",
                        help="Create PSA/FAS comparison")
    parser.add_argument("-c", "--codebase", required=False,
                        dest="codebase", default=None,
                        help="method used for the simulation")
    args = parser.parse_args()

    return args

#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------

ARGS = parse_arguments()

TOP_INPUT_DIR = ARGS.input_dir
if not os.path.isdir(TOP_INPUT_DIR):
    print("Invalid input directory!")
    sys.exit(-1)
if not "Sims" in os.listdir(TOP_INPUT_DIR):
    print("Please provide the top-level simulation directory!\n"
          "This is the directory given to the cluster script")
    sys.exit(-1)
INPUT_OUTDATA = os.path.join(TOP_INPUT_DIR, "Sims", "outdata")

compare_runs(INPUT_OUTDATA, ARGS)
