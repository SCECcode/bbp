#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Broadband Platform script to generate GMPE comparison plot
against ground observations for stations in the station list.
$Id: gmpe_comparison.py 1719 2016-08-18 21:44:13Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

# Import Broadband modules
import bband_utils
import gmpe_config
from PlotGOF import PlotGoF
from install_cfg import InstallCfg
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

class GMPEComparison(object):
    """
    This class implements the methods used to generate the GMPE data
    for a list of stations and then compare the results against ground
    observations files located in a given directory.
    """

    def __init__(self, i_r_stations, i_r_src_file, i_comp_label,
                 i_gmpe_group_name, sim_id=0):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.r_src_file = i_r_src_file
        self.comp_label = i_comp_label
        self.gmpe_group_name = i_gmpe_group_name
        self.src_keys = None

    def read_gmpe(self, input_file):
        """
        Reads the GMPE input_file and returns the periods along with the data
        """
        # Start empty
        gmpe_data = []
        gmpe_models = None

        gmpefile = open(input_file, 'r')
        for line in gmpefile:
            line = line.strip()
            if line.startswith('#period'):
                # This line contains the models we want to use
                gmpe_models = line.split(None, 1)[1]
                gmpe_models = gmpe_models.split()
            if line.startswith('#'):
                # Skip other comments
                continue
            values = line.split()
            values = [float(value) for value in values]
            # Extract period and rotd50 value
            period = values[0]
            medians = values[1:]
            gmpe_data.append((period, medians))
        gmpefile.close()

        # Make sure we parsed the line with the model names
        if gmpe_models is None:
            raise bband_utils.ProcessingError("Cannot find GMPE models in %s" %
                                              (input_file))

        # Return the station median data
        return gmpe_data, gmpe_models

    def read_rotd50(self, input_file):
        """
        Reads the RotD50 input_file and returns the periods with rotd50 data
        """
        # Start empty
        periods = []
        rd50 = []

        rd50file = open(input_file, 'r')
        for line in rd50file:
            line = line.strip()
            if line.startswith('#'):
                # Skip comments
                continue
            values = line.split()
            if len(values) != 4:
                # We are looking for 4 items, skip this line
                continue
            # Convert to floats
            values = [float(value) for value in values]
            # Extract period and rotd50 value
            periods.append(values[0])
            rd50.append(values[3])
        rd50file.close()

        return (periods, rd50)

    def calculate_residuals(self, station, gmpe_model, gmpe_data,
                            obs_periods, obs_data, resid_file,
                            print_headers):
        """
        This function calculates the residuals for the gmpe data
        versus the obs_data, and outputs the results to the resid_file
        """
        # Get gmpe periods
        gmpe_periods = [points[0] for points in gmpe_data]
        # Find common set
        period_set = sorted(list(set(gmpe_periods).intersection(obs_periods)))
        # Create new index arrays for observations and gmpes
        gmpe_items = []
        obs_items = []
        for period in period_set:
            gmpe_items.append(gmpe_periods.index(period))
            obs_items.append(obs_periods.index(period))
        # Get gmpe data array
        gmpe_group = gmpe_config.GMPES[self.gmpe_group_name]
        index = gmpe_group["models"].index(gmpe_model)
        res1 = [points[1][index] for points in gmpe_data]

        # Update gmpe_data, and obs_data arrays
        gmpe_points = []
        obs_points = []
        for item1, item2 in zip(gmpe_items, obs_items):
            gmpe_points.append(res1[item1])
            obs_points.append(obs_data[item2])

        # Calculate residuals
        for idx in range(0, len(obs_points)):
            if gmpe_points[idx] != 0.0:
                gmpe_points[idx] = math.log(obs_points[idx]/gmpe_points[idx])
            else:
                gmpe_points[idx] = -99

        # Now, output to file
        if print_headers:
            outf = open(resid_file, 'w')
            outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                       ("EQ", "Mag", "stat", "lon", "lat", "stat_seq_no",
                        "Vs30", "close_dist", "Xcos", "Ycos", "T_min",
                        "T_max", "comp"))
            for period in period_set:
                outf.write("\t%.5e" % (period))
            outf.write("\n")
        else:
            outf = open(resid_file, 'a')

        # Calculate Rrup
        origin = (self.src_keys['lon_top_center'],
                  self.src_keys['lat_top_center'])
        dims = (self.src_keys['fault_length'], self.src_keys['dlen'],
                self.src_keys['fault_width'], self.src_keys['dwid'],
                self.src_keys['depth_to_top'])
        mech = (self.src_keys['strike'], self.src_keys['dip'],
                self.src_keys['rake'])

        site_geom = [float(station.lon), float(station.lat), 0.0]
        (fault_trace1, up_seis_depth,
         low_seis_depth, ave_dip,
         dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
        _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                         fault_trace1,
                                                         up_seis_depth,
                                                         low_seis_depth,
                                                         ave_dip)

        outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                   (self.comp_label, str(self.src_keys['magnitude']),
                    station.scode, station.lon, station.lat, "-999",
                    station.vs30, rrup, "-999", "-999"))

        if station.high_freq_corner > 0:
            outf.write("\t%.3f" %
                       (1.0 / station.high_freq_corner))
        else:
            outf.write("\t-99999.999")
        if station.low_freq_corner > 0:
            outf.write("\t%.3f" %
                       (1.0 / station.low_freq_corner))
        else:
            outf.write("\t-99999.999")
        outf.write("\t%s" % (gmpe_model.lower()))
        for value in gmpe_points:
            outf.write("\t%.5e" % (value))
        outf.write("\n")
        outf.close()

        return period_set

    def run(self):
        """
        Calculate GMPEs, create bias plot comparisons
        """
        print("GMPE Comparison".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        # Input, tmp, and output directories
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir_gmpe = os.path.join(install.A_OUT_DATA_DIR, str(sim_id),
                                     "gmpe_data_%s" % (sta_base))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))

        self.log = os.path.join(a_logdir, "%d.gmpe_compare.log" % (sim_id))

        #
        # Make sure the output and tmp directories exist
        #
        dirs = [a_tmpdir, a_tmpdir_seis, a_outdir_gmpe, a_outdir, a_logdir]
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
        print_headers = True
        gmpe_models = []
        for site in site_list:
            stat = site.scode
            obs_file = os.path.join(a_tmpdir_seis, "%s.rd50" % (stat))
            gmpe_file = os.path.join(a_outdir_gmpe, "%s-gmpe.ri50" % (stat))
            # Skip station if we don't have observation file
            if not os.access(obs_file, os.R_OK):
                continue
            gmpe_data, gmpe_models[:] = self.read_gmpe(gmpe_file)
            obs_periods, obs_data = self.read_rotd50(obs_file)

            # Loop through the NGA methods
            for gmpe_model in gmpe_models:
                resid_file = os.path.join(a_outdir_gmpe, "%s-%d.resid.txt" %
                                          (gmpe_model.lower(), sim_id))
                period_set = self.calculate_residuals(site, gmpe_model,
                                                      gmpe_data, obs_periods,
                                                      obs_data, resid_file,
                                                      print_headers)
            print_headers = False

        for gmpe_model in gmpe_models:
            # Now call the resid2uncer_varN program to summarize the
            # residuals and create the files needed for the GOF plot
            resid_file = os.path.join(a_outdir_gmpe, "%s-%d.resid.txt" %
                                      (gmpe_model.lower(), sim_id))
            fileroot = os.path.join(a_outdir, "%s-GMPE-%d_r%d-all-rd50-%s" %
                                    (self.comp_label, sim_id,
                                     0, gmpe_model.lower()))
            cmd = ("%s " % os.path.join(install.A_GP_BIN_DIR,
                                        "resid2uncer_varN") +
                   "residfile=%s fileroot=%s " % (resid_file, fileroot) +
                   "comp=%s nstat=%d nper=%d " % (gmpe_model.lower(),
                                                  len(site_list),
                                                  len(period_set)) +
                   "min_cdst=%d >> %s 2>&1" %
                   (0, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Plot GOF plot
        gmpe_group = gmpe_config.GMPES[self.gmpe_group_name]
        gmpe_labels = gmpe_group["labels"]
        plotter = PlotGoF()
        plottitle = "Comparison between GMPEs and %s" % (self.comp_label)
        fileroot = "%s-GMPE-%d_r%d-all-rd50-" % (self.comp_label, sim_id, 0)
        dataroot = ["%s%s" % (fileroot, model.lower()) for model in gmpe_models]
        plotter.multi_plot(plottitle, dataroot, a_outdir,
                           a_outdir, gmpe_labels, len(site_list))

        print("GMPE Comparison Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 5:
        print("Usage: %s " % (PROG_BASE) +
              "station_list src_file label sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = GMPEComparison(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                        sim_id=int(sys.argv[5]))
    ME.run()
