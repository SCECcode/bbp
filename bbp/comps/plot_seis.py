#!/usr/bin/python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module is responsible for plotting velocity and acceleration seismograms.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
import plot_seismograms
from station_list import StationList
from install_cfg import InstallCfg

# Import Pynga and its utilities
import pynga.utils as putils

class PlotSeis(object):
    """
    This module can be used to plot seismograms inside the Broadband
    Platform.
    """

    def __init__(self, i_r_stations, i_r_srcfile,
                 plot_vel, plot_acc, sim_id=0):
        """
        Initialize basic class parameters
        """
        self.r_stations = i_r_stations
        self.plot_vel = plot_vel
        self.plot_acc = plot_acc
        self.sim_id = sim_id

        install = InstallCfg.getInstance()
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(self.sim_id))

        if i_r_srcfile is not None and i_r_srcfile != "":
            i_a_srcfile = os.path.join(a_indir, i_r_srcfile)
            self.src_keys = bband_utils.parse_src_file(i_a_srcfile)
        else:
            self.src_keys = None

    def run(self):
        """
        This function generates velocity and/or acceleration
        seismogram plots, as requested by the user
        """
        print("Plot Seismograms".center(80, '-'))

        if not self.plot_vel and not self.plot_acc:
            # Nothing needs to be plotted
            return
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))

        a_statlist = os.path.join(a_indir, self.r_stations)
        slo = StationList(a_statlist)
        site_list = slo.getStationList()

        # Get fault information, if available
        if self.src_keys is not None:
            origin = (self.src_keys['lon_top_center'],
                      self.src_keys['lat_top_center'])
            dims = (self.src_keys['fault_length'], self.src_keys['dlen'],
                    self.src_keys['fault_width'], self.src_keys['dwid'],
                    self.src_keys['depth_to_top'])
            mech = (self.src_keys['strike'], self.src_keys['dip'],
                    self.src_keys['rake'])

        for site in site_list:
            print("==> Plotting station: %s" % (site.scode))
            # Calculate Rrup
            rrup = None
            if self.src_keys is not None:
                site_geom = [float(site.lon), float(site.lat), 0.0]
                (fault_trace1, up_seis_depth,
                 low_seis_depth, ave_dip,
                 dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
                _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                                 fault_trace1,
                                                                 up_seis_depth,
                                                                 low_seis_depth,
                                                                 ave_dip)

            # Check if we need to plot velocity seismograms
            if self.plot_vel:
                print("===> Plotting velocity...")
                filename = os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                        (sim_id, site.scode))
                outfile = os.path.join(a_outdir, "%d.%s_velocity_seis.png" %
                                       (sim_id, site.scode))
                plot_seismograms.plot_seis(site.scode, filename, sim_id,
                                           'vel', outfile,
                                           rrup=rrup)
            # Check if we need to plot acceleration seismograms
            if self.plot_acc:
                print("===> Plotting acceleration...")
                filename = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                        (sim_id, site.scode))
                outfile = os.path.join(a_outdir, "%d.%s_acceleration_seis.png" %
                                       (sim_id, site.scode))
                plot_seismograms.plot_seis(site.scode, filename, sim_id,
                                           'acc', outfile,
                                           rrup=rrup)

        print("Plot Seismograms Completed".center(80, '-'))

if __name__ == '__main__':
    ME = PlotSeis(sys.argv[1], True, True, sim_id=int(sys.argv[2]))
    ME.run()
