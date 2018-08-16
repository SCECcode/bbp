#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Broadband Platform Version of Rob Graves gen_resid_table.csh and resid2uncer.csh
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
from gp_gof_cfg import GPGofCfg
from install_cfg import InstallCfg
from station_list import StationList
from PlotGOF import PlotGoF
from plot_dist_gof import plot_dist_gof
from plot_map_gof import plot_map_gof

# Import Pynga and its utilities
import pynga.utils as putils

class GPGof(object):
    """
    This class generates GOF plots for the rotd50 data
    """

    def __init__(self, i_r_srcfile, i_r_stations,
                 i_mag, i_comparison_label, cutoff=None,
                 single_component=False, sim_id=0):
        """
        Initialize class instance variables
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_stations = i_r_stations
        self.mag = i_mag
        self.comp_label = i_comparison_label
        self.max_cutoff = cutoff
        self.install = None
        self.config = None
        self.src_keys = None
        if single_component:
            self.single_component = True
        else:
            self.single_component = False

    def summarize_rotd50(self, site_list, a_outdir, a_outdir_gmpe):
        """
        Summarizes all rotd50 data and creates the rotd50 GOF plot
        """
        sim_id = self.sim_id
        install = self.install
        config = self.config

        rd50_residfile = os.path.join(a_outdir, "%s-%d.rd50-resid.txt" %
                                      (self.comp_label, sim_id))
        for comp in config.COMPS_PSA5:
            # Build paths and check lengths
            fileroot = os.path.join(a_outdir, "%s-%d_r%d-%d-rd50-%s" %
                                    (self.comp_label, sim_id, config.MIN_CDST,
                                     self.max_cutoff, comp))
            bband_utils.check_path_lengths([rd50_residfile, fileroot],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/resid2uncer_varN " % (install.A_GP_BIN_DIR) +
                   "residfile=%s fileroot=%s " % (rd50_residfile, fileroot) +
                   "comp=%s nstat=%d nper=63 " % (comp, len(site_list)) +
                   "min_cdst=%d max_cdst=%d >> %s 2>&1" %
                   (config.MIN_CDST, self.max_cutoff, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Plot GOF for psa5/rotd50 data
        if self.single_component:
            plot_mode = 'rd50-single'
        else:
            plot_mode = 'rd50'
        fileroot = '%s-%d_r0-%d-rd50' % (self.comp_label, sim_id, self.max_cutoff)
        plottitle = ("GOF Comparison between %s and simulation %d" %
                     (self.comp_label, sim_id))
        plotter = PlotGoF()
        plotter.plot(plottitle, fileroot, a_outdir, a_outdir,
                     cutoff=self.max_cutoff, mode=plot_mode, colorset='single')

        # Finally, plot the distance and map GOFs
        plot_dist_gof(rd50_residfile, self.comp_label,
                      a_outdir_gmpe, sim_id=self.sim_id)
        plot_map_gof(self.r_srcfile, self.r_stations, rd50_residfile,
                     self.comp_label, sim_id=self.sim_id)

    def run(self):
        """
        This function in the main entry point for this module. It runs
        the gp_gof component.
        """
        print("GP GoF".center(80, '-'))

        # Initialize basic variables
        self.install = InstallCfg.getInstance()
        self.config = GPGofCfg()
        install = self.install
        config = self.config
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.gp_gof.log" %
                                (sim_id))

        # Input, tmp, and output directories
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_seis = os.path.join(install.A_OUT_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir_gmpe = os.path.join(install.A_OUT_DATA_DIR, str(sim_id),
                                     "gmpe_data_%s" % (sta_base))

        # Source file, parse it!
        a_srcfile = os.path.join(install.A_IN_DATA_DIR,
                                 str(sim_id),
                                 self.r_srcfile)
        self.src_keys = bband_utils.parse_src_file(a_srcfile)

        # Station file
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        # List of observed seismogram files
        filelist = os.listdir(a_outdir_seis)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # check cutoff value
        if self.max_cutoff is None:
            self.max_cutoff = config.MAX_CDST

        print_header_rd50 = 1
        # Remove rd50 resid file
        rd50_resid_output = os.path.join(a_outdir, "%s-%d.rd50-resid.txt" %
                                         (self.comp_label, sim_id))
        if os.path.exists(rd50_resid_output):
            os.remove(rd50_resid_output)

        for site in site_list:
            slon = float(site.lon)
            slat = float(site.lat)
            stat = site.scode

            # Now process rd50 files
            expected_rd50_file = os.path.join(a_outdir, "%d.%s.rd50" %
                                              (sim_id, stat))
            if not os.path.exists(expected_rd50_file):
                # just skip it
                print("Skipping rotd50/psa5 for station %s..." % (stat))
                continue

            # See if the rd50 file exist for comparison. If it doesn't
            # exist, skip this station
            rd50_file = None
            if ("%s.rd50" % (stat)) in filelist:
                rd50_file = "%s.rd50" % (stat)
            else:
                # Skip this station
                continue

            # Calculate Rrup
            origin = (self.src_keys['lon_top_center'],
                      self.src_keys['lat_top_center'])
            dims = (self.src_keys['fault_length'], self.src_keys['dlen'],
                    self.src_keys['fault_width'], self.src_keys['dwid'],
                    self.src_keys['depth_to_top'])
            mech = (self.src_keys['strike'], self.src_keys['dip'],
                    self.src_keys['rake'])

            site_geom = [float(site.lon), float(site.lat), 0.0]
            (fault_trace1, up_seis_depth,
             low_seis_depth, ave_dip,
             dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
            _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                             fault_trace1,
                                                             up_seis_depth,
                                                             low_seis_depth,
                                                             ave_dip)

            # Create path names and check if their sizes are within bounds
            datafile1 = os.path.join(a_outdir_seis, rd50_file)
            simfile1 = os.path.join(a_outdir, "%d.%s.rd50" %
                                    (sim_id, stat))
            outfile = os.path.join(a_outdir, "%s-%d.rd50-resid.txt" %
                                   (self.comp_label, self.sim_id))
            bband_utils.check_path_lengths([datafile1, simfile1, outfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/gen_resid_tbl_3comp bbp_format=1 " %
                   (install.A_GP_BIN_DIR) +
                   "datafile1=%s simfile1=%s " % (datafile1, simfile1) +
                   "comp1=psa5n comp2=psa5e comp3=rotd50 " +
                   "eqname=%s mag=%s stat=%s lon=%.4f lat=%.4f " %
                   (self.comp_label, self.mag, stat, slon, slat) +
                   "vs30=%d cd=%.2f " % (site.vs30, rrup) +
                   "flo=%f fhi=%f " % (site.low_freq_corner,
                                       site.high_freq_corner) +
                   "print_header=%d >> %s 2>> %s" %
                   (print_header_rd50, outfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Only need to print header the first time
            if print_header_rd50 == 1:
                print_header_rd50 = 0

        # Finished per station processing, now summarize and plot the data
        if os.path.exists(rd50_resid_output):
            self.summarize_rotd50(site_list, a_outdir, a_outdir_gmpe)

        print("GP GoF Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 8:
        print("Usage: %s " % (PROG_BASE) +
              "source_file station_list magnitude "
              "comp_label cut_off sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = GPGof(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
               cutoff=int(sys.argv[5]),
               single_component=int(sys.argv[6]),
               sim_id=int(sys.argv[7]))
    ME.run()
    sys.exit(0)
