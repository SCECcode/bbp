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

Program to merge a set of BBP scenario simulations into a single run
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import math
import shutil
import optparse

def copy_sim_data(sim_in_dir, sim_out_dir,
                  new_in_dir, new_out_dir,
                  sim, new_sim, sta_base):
    """
    Copy the data from sim_out_dir to new_out_dir, renaming files as
    needed. Also copy the station list from sim_in_dir, and output it
    to the combined sta_base station list in new_in_dir
    """
    # Combined station file
    output_station_list = os.path.join(new_in_dir, sta_base)
    comb_sta = open(output_station_list, 'a')

    station_list = glob.glob("%s/*.stl" % (sim_in_dir))[0]
    station_list = os.path.join(sim_in_dir, station_list)
    sta_file = open(station_list, 'r')
    stations = []
    for line in sta_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("%") or line.startswith("#"):
            continue
        # Write station line to output combined file
        comb_sta.write("%s\n" % (line))
        # Remember station name
        stations.append(line.split()[2])
    # Close input station list
    sta_file.close()
    # Close combined station file
    comb_sta.close()

    # Now loop through all stations
    for station in stations:
        # Copy files
        for ftype in [".acc.bbp", ".vel.bbp", ".rd50",
                      "_acceleration_seis.png", "_velocity_seis.png"]:
            src_file = os.path.join(sim_out_dir,
                                    "%s.%s%s" % (sim, station, ftype))
            dst_file = os.path.join(new_out_dir,
                                    "%s.%s%s" % (new_sim, station, ftype))
            shutil.copy2(src_file, dst_file)
    # Now figure out if we need to copy other files
    if not glob.glob("%s/*.srf" % (new_out_dir)):
        srf_file = glob.glob("%s/*.srf" % sim_out_dir)[0]
        dst_file = os.path.join(new_out_dir,
                                os.path.basename(srf_file))
        shutil.copy2(src_file, dst_file)

def main():
    """
    Parse command-line options
    """
    prog_base = os.path.basename(sys.argv[0])
    usage = "usage: %s [options]" % (prog_base)
    parser = optparse.OptionParser(usage)
    parser.add_option("-d", "--dir", type="string", action="store",
                      dest="simdir",
                      help="Simulation directory")
    parser.add_option("-s", "--sim-id", type="string", action="store",
                      dest="new_sim", help="Simulation id for the merged runs")

    (options, _) = parser.parse_args()

    simdir = options.simdir
    if simdir is None:
        print("Please provide a simulation directory!")
        sys.exit(1)
    simdir = os.path.abspath(simdir)

    new_sim = options.new_sim

    indir = os.path.join(simdir, "Sims", "indata")
    outdir = os.path.join(simdir, "Sims", "outdata")
    newindir = os.path.join(indir, new_sim)
    newoutdir = os.path.join(outdir, new_sim)

    one_sim = os.listdir(indir)[0]
    one_sim = os.path.join(indir, one_sim)
    sta_base = glob.glob("%s/*.stl" % (one_sim))[0]
    sta_base = os.path.basename(sta_base)
    sta_base = "%s.stl" % (sta_base[:sta_base.rfind("-")])

    sims = os.listdir(outdir)
    # Create directories for new simulation
    if not os.path.exists(newindir):
        os.makedirs(newindir)
    if not os.path.exists(newoutdir):
        os.makedirs(newoutdir)

    # Loop through the simulations
    for sim in sims:
        sim_in_dir = os.path.join(indir, sim)
        sim_out_dir = os.path.join(outdir, sim)
        copy_sim_data(sim_in_dir, sim_out_dir,
                      newindir, newoutdir,
                      sim, new_sim, sta_base)
    # All done!

if __name__ == "__main__":
    main()
