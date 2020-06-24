#!/usr/bin/env python
"""
Copyright 2010-2020 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Gen_HTML module for html and index generation
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import time
import shutil

# Import Broadband modules
import bband_utils
import validation_cfg
import velocity_models
from install_cfg import InstallCfg
from station_list import StationList

class GenHTML(object):
    """
    Implement html generation as python component
    """

    def __init__(self, i_r_stations, i_r_src_file, i_vmodel_name,
                 i_val_name, i_method, sim_id=0):
        """
        Initialize class variables
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.r_src_file = i_r_src_file
        self.vmodel_name = i_vmodel_name
        self.val_name = i_val_name
        self.method = i_method
        self.log = None

    def run(self):
        """
        Generate an index file in the outdata directory
        """
        print("GenHTML".center(80, '-'))

        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.genhtml.log" % (sim_id))
        a_statfile = os.path.join(a_indir, self.r_stations)
        a_param_outdir = os.path.join(a_outdir, "param_files")
        a_param_statfile = os.path.join(a_param_outdir, self.r_stations)
        if self.r_src_file is not None and self.r_src_file != "":
            a_src_file = os.path.join(a_indir, self.r_src_file)
            a_param_srcfile = os.path.join(a_param_outdir, self.r_src_file)
            src_props = bband_utils.parse_properties(a_src_file)
            if "seed" in src_props:
                seed = src_props["seed"]
            else:
                seed = "not available"
        else:
            a_src_file = None
            a_param_srcfile = None

        # Make sure tmpdir, outdir exist
        dirs = [a_tmpdir, a_outdir, a_param_outdir]
        bband_utils.mkdirs(dirs, print_cmd=False)

        # Copy station list, srf_file to outdir's param_files directory
        shutil.copy2(a_statfile, a_param_statfile)
        if a_param_srcfile is not None:
            shutil.copy2(a_src_file, a_param_srcfile)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        vel_version = ("%s - %s" % (vel_obj.get_name(), vel_obj.get_version()))

        # Get pointer to validation object, if any
        val_version = None
        if self.val_name:
            val_obj = validation_cfg.VE_EVENTS.get_event_by_name(self.val_name)
            if val_obj is not None:
                val_version = ("%s - %s" % (val_obj.get_print_name(),
                                            val_obj.get_version()))

        #
        # Read and parse the station list with this call
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        index_file = os.path.join(a_outdir, "index-%d.html" % (sim_id))
        idxout = open(index_file, 'w')
        idxout.write("<html>\n")
        idxout.write("<title>Results for simulation %d</title>\n" % (sim_id))
        idxout.write("<body>\n")
        idxout.write("<h2>Simulation Results</h2>\n")
        idxout.write("<table>\n")
        idxout.write("<tr>\n")
        idxout.write("<td>Broadband Version</td>\n")
        idxout.write("<td>%s</td>\n" % (install.VERSION))
        idxout.write("</tr>\n")
        idxout.write("<tr>\n")
        idxout.write("<td>Velocity model version</td>\n")
        idxout.write("<td>%s</td>\n" % (vel_version))
        idxout.write("</tr>\n")
        if val_version:
            idxout.write("<tr>\n")
            idxout.write("<td>Validation package version</td>\n")
            idxout.write("<td>%s</td>\n" % (val_version))
            idxout.write("</tr>\n")
        if install.start_time is not None:
            idxout.write("<tr>\n")
            idxout.write("<td>Simulation Start Time</td>\n")
            idxout.write("<td>%s</td>\n" %
                         (time.strftime("%a %d %b %Y %X %Z",
                                        install.start_time)))
            idxout.write("</tr>\n")
        idxout.write("<tr>\n")
        idxout.write("<td>Simulation End Time</td>\n")
        idxout.write("<td>%s</td>\n" %
                     (time.strftime("%a %d %b %Y %X %Z",
                                    time.localtime())))
        idxout.write("</tr>\n")
        idxout.write("<tr>\n")
        idxout.write("<td>Simulation ID</td>\n")
        idxout.write("<td>%d</td>\n" % (sim_id))
        idxout.write("</tr>\n")
        idxout.write("<tr>\n")
        idxout.write("<td>Simulation Method</td>\n")
        idxout.write("<td>%s</td>\n" % (self.method))
        idxout.write("</tr>\n")
        # Add xml file
        if os.path.exists(os.path.join(a_outdir, "%d.xml" % (sim_id))):
            idxout.write("<tr>\n")
            idxout.write("<td>Sim Spec</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%d.xml" % (sim_id)),
                          "%d.xml" % (sim_id)))
            idxout.write("</tr>\n")
        # Add station list and src_file
        if os.path.exists(os.path.join(a_param_outdir, self.r_stations)):
            idxout.write("<tr>\n")
            idxout.write("<td>Station List</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "param_files", self.r_stations),
                          self.r_stations))
            idxout.write("</tr>\n")
        if a_param_srcfile is not None:
            if os.path.exists(os.path.join(a_param_outdir, self.r_src_file)):
                idxout.write("<tr>\n")
                idxout.write("<td>Source Description</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".",
                                           "param_files",
                                           self.r_src_file),
                              self.r_src_file))
                idxout.write("</tr>\n")
                idxout.write("<tr>\n")
                idxout.write("<td>Random Seed</td>\n")
                idxout.write('<td>%s</td>\n' % (seed))
                idxout.write("</tr>\n")
        # Get bias plots
        dist_lin_plot = glob.glob(os.path.join(a_outdir, "gof-dist-lin*.png"))
        dist_log_plot = glob.glob(os.path.join(a_outdir, "gof-dist-log*.png"))
        rd50plot = glob.glob(os.path.join(a_outdir, "gof*-rd50.png"))
        gmpegofplot = glob.glob(os.path.join(a_outdir, "gof*-GMPE-*.png"))
        mapgofplot = glob.glob(os.path.join(a_outdir, "gof-map-*.png"))
        vs30gofplot = glob.glob(os.path.join(a_outdir, "gof-vs30*.png"))
        if len(gmpegofplot) == 1:
            gmpegofplot = gmpegofplot[0]
        else:
            gmpegofplot = ""
        if len(mapgofplot) == 1:
            mapgofplot = mapgofplot[0]
        else:
            mapgofplot = ""
        if len(vs30gofplot) == 1:
            vs30gofplot = vs30gofplot[0]
        else:
            vs30gofplot = ""
        if len(dist_lin_plot) == 1:
            dist_lin_plot = dist_lin_plot[0]
        else:
            dist_lin_plot = ""
        if len(dist_log_plot) == 1:
            dist_log_plot = dist_log_plot[0]
        else:
            dist_log_plot = ""
        if len(rd50plot) == 1:
            rd50plot = rd50plot[0]
        else:
            if gmpegofplot:
                rd50plot = [plot for plot in rd50plot if plot != gmpegofplot]
            if mapgofplot:
                rd50plot = [plot for plot in rd50plot if plot != mapgofplot]
            if vs30gofplot:
                rd50plot = [plot for plot in rd50plot if plot != vs30gofplot]
            if dist_lin_plot:
                rd50plot = [plot for plot in rd50plot if plot != dist_lin_plot]
            if dist_log_plot:
                rd50plot = [plot for plot in rd50plot if plot != dist_log_plot]
            if len(rd50plot) == 1:
                rd50plot = rd50plot[0]
            else:
                rd50plot = ""
        gmpegofplot = os.path.basename(gmpegofplot)
        mapgofplot = os.path.basename(mapgofplot)
        vs30gofplot = os.path.basename(vs30gofplot)
        rd50plot = os.path.basename(rd50plot)
        dist_lin_plot = os.path.basename(dist_lin_plot)
        dist_log_plot = os.path.basename(dist_log_plot)

        # Add RotD50 bias plot
        if rd50plot:
            idxout.write("<tr>\n")
            idxout.write("<td>RotD50 Bias Plot</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (rd50plot)),
                          "PNG"))
            idxout.write("</tr>\n")
        # Add RotD50 map plot
        if mapgofplot:
            idxout.write("<tr>\n")
            idxout.write("<td>RotD50 Map GOF Plot</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (mapgofplot)),
                          "PNG"))
            idxout.write("</tr>\n")
        # Add the GMPE bias plot
        if gmpegofplot:
            idxout.write("<tr>\n")
            idxout.write("<td>GMPE Comparison Bias Plot</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (gmpegofplot)),
                          "PNG"))
            idxout.write("</tr>\n")
        # Add distance plots
        if dist_lin_plot:
            idxout.write("<tr>\n")
            idxout.write("<td>RotD50 Dist Bias Linear</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (dist_lin_plot)),
                          "PNG"))
            idxout.write("</tr>\n")
        if dist_log_plot:
            idxout.write("<tr>\n")
            idxout.write("<td>RotD50 Dist Bias Log</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (dist_log_plot)),
                          "PNG"))
            idxout.write("</tr>\n")
        if vs30gofplot:
            idxout.write("<tr>\n")
            idxout.write("<td>RotD50 Vs30 GOF Plot</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "%s" % (vs30gofplot)),
                          "PNG"))
            idxout.write("</tr>\n")
        # Add station map
        if os.path.exists(os.path.join(a_outdir, "station_map.png")):
            idxout.write("<tr>\n")
            idxout.write("<td>Station Map</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", "station_map.png"),
                          "PNG"))
            idxout.write("</tr>\n")
        # Now get SRF file and plot
        srfs = glob.glob(os.path.join(a_outdir, "*.srf"))
        if len(srfs) == 1:
            srffile = os.path.basename(srfs[0])
            srfplot = ("%s.png" %
                       (os.path.basename(os.path.splitext(srffile)[0])))
            if not os.path.exists(os.path.join(a_outdir, srfplot)):
                srfplot = ""
        else:
            srffile = ""
            srfplot = ""
        if srffile:
            idxout.write("<tr>\n")
            idxout.write("<td>Rupture file</td>\n")
            idxout.write('<td><a href="%s">%s</a></td>\n' %
                         (os.path.join(".", srffile),
                          "SRF"))
            if srfplot:
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", srfplot),
                              "PNG"))
            idxout.write("</tr>\n")
        idxout.write("</table>\n")
        idxout.write("<p><p>\n")

        for sits in site_list:
            site = sits.scode
            idxout.write("<p>\n")
            idxout.write("<h2>%s</h2>\n" % (site))
            idxout.write("<table>\n")

            # Find all files
            velfile = "%d.%s.vel.bbp" % (sim_id, site)
            velplot = "%d.%s_velocity_seis.png" % (sim_id, site)
            accfile = "%d.%s.acc.bbp" % (sim_id, site)
            accplot = "%d.%s_acceleration_seis.png" % (sim_id, site)
            rd50file = "%d.%s.rd50" % (sim_id, site)
            rd100file = "%d.%s.rd100" % (sim_id, site)
            rd50file_vertical = "%d.%s.rd50.vertical" % (sim_id, site)
            rd100file_vertical = "%d.%s.rd100.vertical" % (sim_id, site)

            # RotD50 Plot
            rd50plot = glob.glob(os.path.join(a_outdir,
                                              "*_%d_%s_rotd50.png" %
                                              (sim_id, site)))
            if len(rd50plot) == 1:
                rd50plot = os.path.basename(rd50plot[0])
            else:
                rd50plot = ""

            # Overlay Plot
            overlayfile = glob.glob(os.path.join(a_outdir,
                                                 "*_%d_%s_overlay.png" %
                                                 (sim_id, site)))
            if len(overlayfile) == 1:
                overlayfile = os.path.basename(overlayfile[0])
            else:
                overlayfile = ""

            # GMPE Plot
            gmpeplot = glob.glob(os.path.join(a_outdir,
                                              "*_%d_%s_gmpe.png" %
                                              (sim_id, site)))
            if len(gmpeplot) == 1:
                gmpeplot = os.path.basename(gmpeplot[0])
            else:
                gmpeplot = ""

            if os.path.exists(os.path.join(a_outdir, velfile)):
                idxout.write("<tr>\n")
                idxout.write("<td>Velocity</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", velfile),
                              "BBP"))
                if os.path.exists(os.path.join(a_outdir, velplot)):
                    idxout.write('<td><a href="%s">%s</a></td>\n' %
                                 (os.path.join(".", velplot),
                                  "PNG"))
                idxout.write("</tr>\n")
            if os.path.exists(os.path.join(a_outdir, accfile)):
                idxout.write("<tr>\n")
                idxout.write("<td>Acceleration</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", accfile),
                              "BBP"))
                if os.path.exists(os.path.join(a_outdir, accplot)):
                    idxout.write('<td><a href="%s">%s</a></td>\n' %
                                 (os.path.join(".", accplot),
                                  "PNG"))
                idxout.write("</tr>\n")
            if os.path.exists(os.path.join(a_outdir, rd50file)):
                idxout.write("<tr>\n")
                idxout.write("<td>RotD50</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", rd50file),
                              "HOR"))
                if os.path.exists(os.path.join(a_outdir, rd50file_vertical)):
                    idxout.write('<td><a href="%s">%s</a></td>\n' %
                                 (os.path.join(".", rd50file_vertical),
                                  "VER"))
                if rd50plot:
                    idxout.write('<td><a href="%s">%s</a></td>\n' %
                                 (os.path.join(".", rd50plot),
                                  "PNG"))
                idxout.write("</tr>\n")
            if os.path.exists(os.path.join(a_outdir, rd100file)):
                idxout.write("<tr>\n")
                idxout.write("<td>RotD100</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", rd100file),
                              "HOR"))
                if os.path.exists(os.path.join(a_outdir, rd100file_vertical)):
                    idxout.write('<td><a href="%s">%s</a></td>\n' %
                                 (os.path.join(".", rd100file_vertical),
                                  "VER"))
                idxout.write("</tr>\n")
            if overlayfile:
                idxout.write("<tr>\n")
                idxout.write("<td>Overlay</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", overlayfile),
                              "PNG"))
                idxout.write("</tr>\n")
            if gmpeplot:
                idxout.write("<tr>\n")
                idxout.write("<td>GMPE Plot</td>\n")
                idxout.write('<td><a href="%s">%s</a></td>\n' %
                             (os.path.join(".", gmpeplot),
                              "PNG"))
                idxout.write("</tr>\n")

            idxout.write("</table>\n")

        idxout.write("</body>\n")
        idxout.write("</html>\n")
        idxout.close()

        print("==> Wrote file: %s" % (index_file))
        print("GenHTML Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = GenHTML(sys.argv[1], sys.argv[2], sys.argv[3],
                 sys.argv[4], sys.argv[5],
                 sim_id=int(sys.argv[6]))
    ME.run()
