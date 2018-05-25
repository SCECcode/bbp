#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Broadband Platform script to generate GMPE data for stations in the
station list.
$Id: calculate_gmpe.py 1719 2016-08-18 21:44:13Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
import gmpe_config
from install_cfg import InstallCfg
from station_list import StationList

# Import PyNGA modules
import pynga.utils as putils

class CalculateGMPE(object):
    """
    This class implements the methods used to generate the GMPE data
    for a list of stations.
    """

    def __init__(self, i_r_stations, i_r_src_file,
                 i_data_corrected, i_gmpe_group_name,
                 sim_id=0):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.r_src_file = i_r_src_file
        self.data_corrected = i_data_corrected
        self.gmpe_group_name = i_gmpe_group_name
        self.src_keys = None

    def calculate_gmpe(self, station, output_file):
        """
        This function calculates the gmpe for a station and writes the
        output in the output_file
        """
        gmpe_group = gmpe_config.GMPES[self.gmpe_group_name]
        src_keys = self.src_keys
        origin = (src_keys['lon_top_center'], src_keys['lat_top_center'])
        dims = (src_keys['fault_length'], src_keys['dlen'],
                src_keys['fault_width'], src_keys['dwid'],
                src_keys['depth_to_top'])
        mech = (src_keys['strike'], src_keys['dip'], src_keys['rake'])

        # Station location
        site_geom = [float(station.lon), float(station.lat), 0.0]
        (fault_trace1, upper_seis_depth,
         lower_seis_depth, ave_dip,
         dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
        rjb, rrup, rx = putils.DistanceToSimpleFaultSurface(site_geom,
                                                            fault_trace1,
                                                            upper_seis_depth,
                                                            lower_seis_depth,
                                                            ave_dip)
        # If using corrected ground observation data, calcualte GMPEs
        # with a Vs30 of 863 m/s
        if self.data_corrected:
            vs30 = 863
        else:
            vs30 = station.vs30
        z10 = None # Let PyNGA calculate it
        z25 = None # Let PyNGA calculate it

        # Compute PSA for this station
        station_median = []
        for period in gmpe_group["periods"]:
            period_medians = []
            for nga_model in gmpe_group["models"]:
                median = gmpe_config.calculate_gmpe(self.gmpe_group_name,
                                                    nga_model,
                                                    src_keys['magnitude'],
                                                    rjb, vs30,
                                                    period,
                                                    rake=src_keys['rake'],
                                                    dip=src_keys['dip'],
                                                    W=src_keys['fault_width'],
                                                    Ztor=src_keys['depth_to_top'],
                                                    Rrup=rrup, Rx=rx,
                                                    Z10=z10, Z25=z25)
                period_medians.append(median)
            station_median.append((period, period_medians))

        # Create label
        file_label = ""
        for nga_model in gmpe_group["models"]:
            file_label = "%s %s" % (file_label, nga_model)
        # Output data to file
        outfile = open(output_file, 'w')
        outfile.write("#station: %s\n" % (station.scode))
        outfile.write("#period%s\n" % (file_label))
        for item in station_median:
            period = item[0]
            vals = item[1]
            out_str = "%.4f" % (period)
            for method in vals:
                out_str = out_str + "\t%.6f" % (method)
            outfile.write("%s\n" % (out_str))
        outfile.close()

        # Return list
        return station_median

    def run(self):
        """
        Calculate GMPEs, create bias plot comparisons
        """
        print("Calculate GMPE".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        # Input, tmp, and output directories
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_gmpe = os.path.join(a_outdir,
                                     "gmpe_data_%s" % (sta_base))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))

        self.log = os.path.join(a_logdir, "%d.gmpe_compare.log" % (sim_id))

        #
        # Make sure the output and tmp directories exist
        #
        dirs = [a_outdir_gmpe, a_outdir, a_logdir]
        bband_utils.mkdirs(dirs, print_cmd=False)

        # Source file, parse it!
        a_srcfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_src_file)
        self.src_keys = bband_utils.parse_src_file(a_srcfile)

        # Station file
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Go through each station, and print comparison headers for
        # the first station we process
        for site in site_list:
            stat = site.scode
            print("==> Calculating GMPE for station: %s" % (stat))
            output_file = os.path.join(a_outdir_gmpe, "%s-gmpe.ri50" % (stat))
            self.calculate_gmpe(site, output_file)

        # All done
        print("Calculate GMPE Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 7:
        print("Usage: %s " % (PROG_BASE) +
              "station_list src_file label corr_flag "
              "gmpe_group sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = CalculateGMPE(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                       sys.argv[5], sim_id=int(sys.argv[6]))
    ME.run()
