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

Utility classes for the SCEC Broadband Platform
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
#import cc
#import ggp
import bband_utils
from station_list import StationList
from geobb_srf import GeoBBSRF

def stationlist2filelist(slo):
    """
    Input a stationlist object (slo). From this extract the station codes and
    create a list of file names that represent long period data files
    """
    site_list = slo.getStationList()
    files = []
    for sta in site_list:
        fname = "%s" % (sta.scode)
        files.append(fname)

    return files

def gp_subset(in_file1, in_format1, in_file2, outfile):
    """
    Takes two input stat files in a given format and outputs
    an intersection file in GP format
    """
    #Get a station list from in_file1
    if in_format1 == 'GP' or in_format1 == 'UCSB':
        stat_list = StationList(in_file1).getStationList()
        stat_names = [stat.scode for stat in stat_list]
    elif in_format1 == 'SDSU':
        stat_file_fp = open(in_file1, 'r')
        data = stat_file_fp.readlines()
        stat_file_fp.close()
        for i in range(0, len(data)):
            pieces = data[i].split()
            if len(pieces) > 1:
                if pieces[1] == 'X':
                    break
        stat_names = []
        for j in range(i + 1, len(data)):
            pieces = data[j].split()
            stat_names.append(pieces[2])
    else:
        raise bband_utils.ParameterError("Format %s is not supported." %
                                         (in_format1))

    # Use the station list to subset in_file2
    intersect_list = []
    stat2_list = StationList(in_file2).getStationList()
    for stat2 in stat2_list:
        for stat1 in stat_names:
            if stat2.scode == stat1:
                intersect_list.append(stat2)
    StationList.build(intersect_list, outfile)

def sdsu_subset(working_dir, stat_list_obj, val_obj, sdsu_out_list,
                vsrealfile="VSREAL", vssynthfile="VSSYNTH"):
    """
    This function selects stations from the SDSU stations list and
    trims the ones that are not also present in the from
    stat_list_obj. It creates three output files in the working_dir
    directory, one with the trimmed station list, and one for each of
    vsreal and vssynth.
    """

    # Assemble pathnames
    a_full_sdsu_stations = val_obj.get_input('SDSU', 'stations')
    a_full_vsreal = val_obj.get_input('SDSU', 'vsreal')
    a_full_vssynth = val_obj.get_input('SDSU', 'vssynth')

    # Create output files
    stl_out = open(sdsu_out_list, 'w')
    vsreal_out = open("%s/%s" % (working_dir, vsrealfile), 'w')
    vssynth_out = open("%s/%s" % (working_dir, vssynthfile), 'w')

    # Read three files (stations, vsreal, vssynth)
    fp = open(a_full_sdsu_stations, 'r')
    full_station_data = fp.readlines()
    fp.close()

    fp = open(a_full_vsreal, 'r')
    full_vsreal_data = fp.readlines()
    fp.close()

    fp = open(a_full_vssynth, 'r')
    full_vssynth_data = fp.readlines()
    fp.close()

    # Initialize dictionaries
    stat_line = dict()
    vs_reals = dict()
    vs_synths = dict()

    # Skip the header
    for i in range(0, len(full_station_data)):
        pieces = full_station_data[i].split()
        stl_out.write(full_station_data[i])
        if len(pieces) > 1:
            if pieces[1] == 'X':
                break

    # Now merge the data from the station list, vsreal, and vssynth
    for j in range(i + 1, len(full_station_data)):
        pieces = full_station_data[j].split()
        #print full_station_data[j]
        stat_line[pieces[2]] = full_station_data[j]
        vs_reals[pieces[2]] = full_vsreal_data[j - i - 1]
        vs_synths[pieces[2]] = full_vssynth_data[j - i - 1]

    # Now select only the stations present also in the stat_list_obj
    for stat in stat_list_obj.getStationList():
        #print stat.scode
        if stat.scode in stat_line:
            stl_out.write(stat_line[stat.scode])
            vsreal_out.write(vs_reals[stat.scode])
            vssynth_out.write(vs_synths[stat.scode])

    # Done, flush and close everything
    stl_out.flush()
    stl_out.close()
    vsreal_out.flush()
    vsreal_out.close()
    vssynth_out.flush()
    vssynth_out.close()

def bbp2sdsu_statlist(working_dir, slo, sdsu_stalist,
                      srf_file, xyz_srf_file, extended_fault,
                      tmpfile="station.coords"):
    """
    Takes BBP station list object and writes SDSU station list
    """

    # Convert station coordinates to cartesian coordinates
    a_srf_file = os.path.join(working_dir, srf_file)
    a_xyz_srf_file = os.path.join(working_dir, xyz_srf_file)
    out_coords = os.path.join(working_dir, tmpfile)

    gbb = GeoBBSRF()
    gbb.run(slo, out_coords, extended_fault,
            a_srf_file, a_xyz_srf_file, 'y')
    if not os.path.exists(out_coords):
        raise bband_utils.ProcessingError("Error converting station coordinates"
                                          " to SDSU format, exiting.")
    coords_fp = open(out_coords, 'r')
    coords_data = coords_fp.readlines()
    coords_fp.close()

    # Write station file
    stat_list = slo.getStationList()
    stalist_fp = open(sdsu_stalist, 'w')
    for i in range(0, len(stat_list)):
        pieces = coords_data[i].split()
        stalist_fp.write("%f  %f  %s  -1  -1  -1\n" % (float(pieces[0]),
                                                       float(pieces[1]),
                                                       stat_list[i].scode))
    stalist_fp.flush()
    stalist_fp.close()

    # All done!
    return "%s.param" % out_coords

def bbp2sdsu_vsfiles(working_dir, slo, vsrealfile, vssynthfile):
    """
    Produces VSREAL from station list and VSSYNTH (always 865)
    """
    stat_list = slo.getStationList()

    # Write VSREAL and VSSYNTH files in working_dir
    if vsrealfile is not None:
        vsreal = os.path.join(working_dir, vsrealfile)
        vsreal_fp = open(vsreal, 'w')
        for site in stat_list:
            vsreal_fp.write("  %d\n" % site.vs30)
        vsreal_fp.flush()
        vsreal_fp.close()

    if vssynthfile is not None:
        vssynth = os.path.join(working_dir, vssynthfile)
        vssynth_fp = open(vssynth, 'w')
        for site in stat_list:
            vssynth_fp.write("   865\n")
        vssynth_fp.flush()
        vssynth_fp.close()

def gp2uc_stalist(slo, uc_stalist, uc_vs30):
    """
    Pass in a station list object (slo) from the standard GP station list.
    All comments have been removed and each line in the list is a station.
    Convert this into a UCSB format station list containing lat/lon values.
    Values in list are strings at this time.
    """
    sl = slo.getStationList()
    f = open(uc_stalist, "w")
    vs = open(uc_vs30, "w")
    num = len(sl)
    my_str = "%d\n" % (num)
    f.write(my_str)
    for i in sl:
        my_str = "%s %s %s\n" % (i.lon, i.lat, i.scode)
        f.write(my_str)
        my_str = "\t%s\t%s\n" % (i.scode, i.vs30)
        vs.write(my_str)
    f.flush()
    f.close()
    vs.flush()
    vs.close()
    return

def sdsu2uc_subset(sdsu_stat_file, uc_stalist_in, uc_vs30_in,
                   uc_stalist_out, uc_vs30_out):
    # Read input file
    stat_file_fp = open(sdsu_stat_file, "r")
    data = stat_file_fp.readlines()
    stat_file_fp.close()
    for i in range(0, len(data)):
        pieces = data[i].split()
        if len(pieces) > 1:
            if pieces[1] == "X":
                break
    stat_names = []
    for j in range(i + 1, len(data)):
        pieces = data[j].split()
        stat_names.append(pieces[2])

    new_list = []
    slo = StationList(uc_stalist_in).getStationList()
    for stat in stat_names:
        for entry in slo:
            if stat == entry.scode:
                new_list.append(entry)
    StationList.build(new_list, uc_stalist_out)

    fp = open(uc_vs30_in, 'r')
    vs30_dict = dict()
    for line in fp.readlines():
        pieces = line.split()
        vs30_dict[pieces[0]] = pieces[1]
    fp.close()

    fp = open(uc_vs30_out, 'w')
    for stat in stat_names:
        if stat in vs30_dict:
            fp.write("\t%s\t%s\n" % (stat, vs30_dict[stat]))
    fp.flush()
    fp.close()

if __name__ == "__main__":
    print("Testing: %s" % (sys.argv[0]))
