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

Program to create the RotD100 comparison
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import install_cfg
import bbp_formatter
from correct_psa import CorrectPSA
from station_list import StationList
from PlotGOF import PlotGoF

# Import Pynga and its utilities
import pynga.utils as putils

SUPPORTED_OBS_FORMATS = ["acc_bbp", "acc_peer"]

def do_rotd100(workdir, peer_input_e_file, peer_input_n_file,
               output_rotd100_file, logfile):
    """
    This function runs the rotd100 command inside workdir, using
    the inputs and outputs specified
    """
    install = install_cfg.InstallCfg.getInstance()

    # Make sure we don't have absolute path names
    peer_input_e_file = os.path.basename(peer_input_e_file)
    peer_input_n_file = os.path.basename(peer_input_n_file)
    output_rotd100_file = os.path.basename(output_rotd100_file)

    # Save cwd, change back to it at the end
    old_cwd = os.getcwd()
    os.chdir(workdir)

    # Make sure we remove the output files first or Fortran will
    # complain if they already exist
    try:
        os.unlink(output_rotd100_file)
    except OSError:
        pass

    #
    # write config file for rotd100 program
    rd100_conf = open("rotd100_inp.cfg", 'w')
    # This flag indicates inputs acceleration
    rd100_conf.write("2 interp flag\n")
    # This flag indicate we are processing two input files
    rd100_conf.write("1 Npairs\n")
    # Number of headers in the file
    rd100_conf.write("6 Nhead\n")
    rd100_conf.write("%s\n" % peer_input_e_file)
    rd100_conf.write("%s\n" % peer_input_n_file)
    rd100_conf.write("%s\n" % output_rotd100_file)
    # Close file
    rd100_conf.close()

    progstring = ("%s >> %s 2>&1" %
                  (os.path.join(install.A_UCB_BIN_DIR,
                                "rotd100"), logfile))
    bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

    # Restore working directory
    os.chdir(old_cwd)

def do_split_rotd50(a_tmpdir, in_rotd100_base, out_rotd50_base, logfile):
    """
    Create a RotD50 file out of the RotD100 file for modules expecting this
    """
    # Open input and output files
    in_file = open(os.path.join(a_tmpdir, in_rotd100_base), 'r')
    out_file = open(os.path.join(a_tmpdir, out_rotd50_base), 'w')

    for line in in_file:
        line = line.strip()
        # Skip comments
        if line.startswith('#'):
            out_file.write("%s\n" % (line))
            continue
        pieces = line.split()
        pieces = pieces[0:4]
        out_file.write("  %s\n" % (" ".join(pieces)))

    # Close everything
    in_file.close()
    out_file.close()

class RotD100(object):
    """
    BBP module implementation of rotd100 provided by UCB.
    Rotd100 inputs seismograms and outputs response spectra
    """

    def __init__(self, i_r_stations, i_r_srcfile=None,
                 i_a_obsdir=None, i_obs_format=None, i_obs_corr=None,
                 i_comparison_label=None, cutoff=None, sim_id=0):
        """
        Initializes class variables
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_stations = i_r_stations
        self.comp_label = i_comparison_label
        self.max_cutoff = cutoff
        self.a_obsdir = i_a_obsdir
        self.obs_format = i_obs_format
        self.obs_corrections = i_obs_corr
        self.src_keys = None

        # Make observed seismograms are in a format we can handle
        if i_obs_format is not None and i_obs_format not in SUPPORTED_OBS_FORMATS:
            raise bband_utils.ParameterError("Format %s for " %
                                             (self.obs_format) +
                                             "observed seismograms "
                                             "not supported")

    def calculate_simulated(self, a_statfile, a_tmpdir, a_outdir, a_dstdir):
        """
        This function calculates the RotD100/RotD50 values for the
        computed seismograms
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for site in site_list:
            stat = site.scode
            print("==> Calculating simulation RotD100 for station: %s" %
                  (stat))
            # Since we have velocity files, we need to differentiate
            # to get to acceleration

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, stat))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, stat))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, stat))
            bbpfile = os.path.join(a_outdir,
                                   "%d.%s.vel.bbp" % (sim_id, stat))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                   "wcc2bbp=0 nsfile=%s ewfile=%s udfile=%s < %s >> %s 2>&1" %
                   (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            for component in ["090", "000", "ver"]:
                # Differentiate to get from velocity to accl needed by rotd100
                # Create path names and check if their sizes are within bounds
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.%s" %
                                      (sim_id, stat, component))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.acc.%s" %
                                       (sim_id, stat, component))

                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s diff=1 " %
                       (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                       "filein=%s fileout=%s >> %s 2>&1" %
                       (filein, fileout, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                # Check file length
                bband_utils.check_path_lengths(["%s/%d.%s.acc.%s" %
                                                (a_tmpdir, sim_id,
                                                 stat, component)],
                                               bband_utils.GP_MAX_FILENAME)

            # Now we need to convert them back to bbp
            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.000" % (sim_id, stat))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.090" % (sim_id, stat))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.ver" % (sim_id, stat))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.acc.bbp" % (sim_id, stat))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Now we need to convert to peer format
            out_n_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_n.acc" % (sim_id, stat))
            out_e_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_e.acc" % (sim_id, stat))
            out_z_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_z.acc" % (sim_id, stat))
            bbp_formatter.bbp2peer(bbpfile, out_n_acc, out_e_acc, out_z_acc)

            # Let's have rotD100 create these output files
            out_rotd50_base = "%d.%s.rd50" % (sim_id, stat)
            out_rotd100_base = "%d.%s.rd100" % (sim_id, stat)
            out_rotd50v_base = "%d.%s.rd50.vertical" % (sim_id, stat)
            out_rotd100v_base = "%d.%s.rd100.vertical" % (sim_id, stat)

            # Run the rotD100 program twice (horizontals and vertical)
            do_rotd100(a_tmpdir, out_e_acc, out_n_acc,
                       out_rotd100_base, self.log)
            # Run the rotD100 program twice (horizontals and vertical)
            do_rotd100(a_tmpdir, out_z_acc, out_z_acc,
                       out_rotd100v_base, self.log)
            # Create rotd50 files as well
            do_split_rotd50(a_tmpdir, out_rotd100_base,
                            out_rotd50_base, self.log)
            do_split_rotd50(a_tmpdir, out_rotd100v_base,
                            out_rotd50v_base, self.log)

            # Copy horizontals
            cmd = "cp %s %s" % (os.path.join(a_tmpdir, out_rotd100_base),
                                os.path.join(a_dstdir, out_rotd100_base))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
            cmd = "cp %s %s" % (os.path.join(a_tmpdir, out_rotd50_base),
                                os.path.join(a_dstdir, out_rotd50_base))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
            # Now copy verticals
            cmd = "cp %s %s" % (os.path.join(a_tmpdir, out_rotd100v_base),
                                os.path.join(a_dstdir, out_rotd100v_base))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
            cmd = "cp %s %s" % (os.path.join(a_tmpdir, out_rotd50v_base),
                                os.path.join(a_dstdir, out_rotd50v_base))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

    def calculate_observations(self, a_indir, a_statfile, a_tmpdir_seis, a_dstdir):
        """
        This function calculates RotD100/RotD50 for the observation
        seismograms. It corrects the observations using the user-provided
        correction coefficients.
        """
        sim_id = self.sim_id
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Inialize the CorrectPSA module
        if self.obs_corrections:
            corr_psa = CorrectPSA(self.r_stations,
                                  "rd100",
                                  os.path.join(a_indir,
                                               self.obs_corrections),
                                  a_tmpdir_seis, sim_id)
        else:
            corr_psa = None

        # List of observed seismogram files
        filelist = os.listdir(self.a_obsdir)

        # Go through each station
        for site in site_list:
            stat = site.scode
            print("==> Calculating observations RotD100 for station: %s" %
                  (stat))
            # Check if we have the corresponding calculated seismogram
            expected_calculated_file = os.path.join(a_dstdir,
                                                    "%d.%s.rd100" %
                                                    (sim_id, stat))
            if not os.path.exists(expected_calculated_file):
                # Just skip it
                print("Couldn't find file: %s" %
                      (expected_calculated_file) +
                      "This is not necessarily an error, as you may have " +
                      "run with a subset of a stations. Continuing " +
                      "with available stations.")
                continue

            # Ok, we have a simulated seismogram for this station,
            # let's look for the observed file
            r_e_peer_file = None
            r_n_peer_file = None
            r_z_peer_file = None
            r_bbp_file = "%s.bbp" % (stat)
            # Do different things depending on the format of the
            # observed seismograms
            if self.obs_format == "acc_bbp":
                # We need to look for the bbp file
                if r_bbp_file not in filelist:
                    # No bbp file for this station
                    continue
                print(r_bbp_file)
                # Copy bbp file to the tmp seismogram directory
                a_src_bbp_file = os.path.join(self.a_obsdir, r_bbp_file)
                a_dst_bbp_file = os.path.join(a_tmpdir_seis, r_bbp_file)
                shutil.copy2(a_src_bbp_file, a_dst_bbp_file)
                # Now we need to create the peer files to process with rotd50
                r_e_peer_file = os.path.join(a_tmpdir_seis, "%s_E.acc" % (stat))
                r_n_peer_file = os.path.join(a_tmpdir_seis, "%s_N.acc" % (stat))
                r_z_peer_file = os.path.join(a_tmpdir_seis, "%s_Z.acc" % (stat))
                bbp_formatter.bbp2peer(a_dst_bbp_file,
                                       r_n_peer_file,
                                       r_e_peer_file,
                                       r_z_peer_file)
            elif self.obs_format == "acc_peer":
                # Look for the E, N, and Z files
                for my_file in filelist:
                    if my_file.endswith("%s_E.acc" % (stat)):
                        r_e_peer_file = my_file
                        if (r_n_peer_file is not None and
                            r_z_peer_file is not None):
                            break
                    elif my_file.endswith("%s_N.acc" % (stat)):
                        r_n_peer_file = my_file
                        if (r_e_peer_file is not None and
                            r_z_peer_file is not None):
                            break
                    elif my_file.endswith("%s_Z.acc" % (stat)):
                        r_z_peer_file = my_file
                        if (r_e_peer_file is not None and
                            r_n_peer_file is not None):
                            break
                if ((r_e_peer_file is None) or
                    (r_n_peer_file is None) or
                    (r_z_peer_file is None)):
                    # Couldn't find all 3 files
                    continue
                #print(r_e_peer_file, r_n_peer_file, r_z_peer_file)
                # Copy all three files to the tmp seismogram directory
                for eachfile in (r_e_peer_file, r_n_peer_file, r_z_peer_file):
                    a_src_peer_file = os.path.join(self.a_obsdir, eachfile)
                    a_dst_peer_file = os.path.join(a_tmpdir_seis, eachfile)
                    shutil.copy2(a_src_peer_file, a_dst_peer_file)

                # Now we need to convert them into bbp format
                bbp_formatter.peer2bbp(os.path.join(a_tmpdir_seis,
                                                    r_n_peer_file),
                                       os.path.join(a_tmpdir_seis,
                                                    r_e_peer_file),
                                       os.path.join(a_tmpdir_seis,
                                                    r_z_peer_file),
                                       os.path.join(a_tmpdir_seis,
                                                    r_bbp_file))
            else:
                raise bband_utils.ParameterError("Format %s for " %
                                                 (self.obs_format) +
                                                 "observed seismograms "
                                                 "not supported")

            # Run RotD100 on this file
            if corr_psa is not None:
                # First calculate rd100/50 and psa5 files
                do_rotd100(a_tmpdir_seis, r_e_peer_file, r_n_peer_file,
                           "%s-orig.rd100" % (stat), self.log)

                # Now we need to correct the RotD100/RotD50 output
                # using the user-supplied correction factors
                corr_psa.correct_station(stat, "rd100")
            else:
                # Use final names for output files
                do_rotd100(a_tmpdir_seis, r_e_peer_file, r_n_peer_file,
                           "%s.rd100" % (stat), self.log)
            shutil.copy2(os.path.join(a_tmpdir_seis, "%s.rd100" % (stat)),
                         os.path.join(a_dstdir, "%s.rd100" % (stat)))

    def trim_rd100_file(self, a_input_file, a_output_file):
        """
        Trims columns 2 and 3 of the input file so that the output file
        contains only the following: period, rotd50, rotd100, ratio
        """
        # Open input and output files
        in_file = open(a_input_file, 'r')
        out_file = open(a_output_file, 'w')

        for line in in_file:
            line = line.strip()
            # Skip comments
            if line.startswith('#'):
                out_file.write("%s\n" % (line))
                continue
            pieces = line.split()
            pieces = pieces[0:1] + pieces[3:]
            out_file.write("  %s\n" % (" ".join(pieces)))

        # Close everything
        in_file.close()
        out_file.close()

    def calculate_ratio(self, a_filename):
        """
        This function reads a file and calculates the RotD100/RotD50
        ratio, adding an extra column at the end of the file.
        """
        basedir = os.path.dirname(a_filename)
        tmp_filename = os.path.join(basedir, "rd100.tmp")

        # Open input and temp file for writing extra column
        in_file = open(a_filename, 'r')
        out_file = open(tmp_filename, 'w')

        # Add extra column to temp file
        for line in in_file:
            line = line.strip()
            # Write comments
            if line.startswith('#'):
                out_file.write("%s\n" % (line))
                continue
            pieces = [float(piece) for piece in line.split()]
            ratio = pieces[-1]/pieces[-2]
            out_file.write("  %s %10.5E\n" % (line, ratio))

        # Close files
        in_file.close()
        out_file.close()

        # Move file
        shutil.move(tmp_filename, a_filename)

    def calculate_ratios(self, a_statfile, a_dstdir):
        """
        This function adds an extra column to the rd100 files
        containing the RotD100/RotD50 ratio. It does this for both
        observations and simulated data files.
        """
        sim_id = self.sim_id
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Loop through all stations
        for site in site_list:
            stat = site.scode

            # Simulated data file
            basename = "%d.%s.rd100" % (sim_id, stat)
            filename = os.path.join(a_dstdir, basename)
            self.calculate_ratio(filename)

            # Observation data file
            basename = "%s.rd100" % (stat)
            filename = os.path.join(a_dstdir, basename)
            self.calculate_ratio(filename)

    def calculate_residuals(self, a_statfile, a_dstdir):
        """
        This function calculates the residuals comparing observations and
        calculated data.
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        print_header = 1
        rd100_resid_output = os.path.join(a_dstdir, "%s-%d-resid-rd100.txt" %
                                          (self.comp_label, sim_id))
        # If output file exists, delete it
        if os.path.exists(rd100_resid_output):
            os.remove(rd100_resid_output)

        # Filenames for tmp files, check if filename size within bounds
        obsfile = os.path.join(a_dstdir, "rd100_obs.txt")
        simfile = os.path.join(a_dstdir, "rd100_sim.txt")
        bband_utils.check_path_lengths([obsfile, simfile,
                                        rd100_resid_output],
                                       bband_utils.GP_MAX_FILENAME)

        # Loop through all stations
        for site in site_list:
            slon = float(site.lon)
            slat = float(site.lat)
            stat = site.scode

            # Trim files
            self.trim_rd100_file(os.path.join(a_dstdir,
                                              "%d.%s.rd100" % (sim_id, stat)),
                                 simfile)
            self.trim_rd100_file(os.path.join(a_dstdir,
                                              "%s.rd100" % (stat)),
                                 obsfile)

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

            cmd = ("%s bbp_format=1 " %
                   (os.path.join(install.A_GP_BIN_DIR,
                                 "gen_resid_tbl_3comp")) +
                   "datafile1=%s simfile1=%s " % (obsfile, simfile) +
                   "comp1=rotd50 comp2=rotd100 comp3=ratio " +
                   "eqname=%s mag=%s stat=%s lon=%.4f lat=%.4f " %
                   (self.comp_label,
                    self.src_keys['magnitude'], stat, slon, slat) +
                   "vs30=%d cd=%.2f " % (site.vs30, rrup) +
                   "flo=%f fhi=%f " % (site.low_freq_corner,
                                       site.high_freq_corner) +
                   "print_header=%d >> %s 2>> %s" %
                   (print_header, rd100_resid_output, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Only need to print header the first time
            if print_header == 1:
                print_header = 0

        # Remove temp files
        try:
            os.remove(obsfile)
            os.remove(simfile)
        except:
            pass

    def generate_plot(self, a_statfile, a_dstdir):
        """
        This function generates the bias plot with ratio of maximum to
        median response across orientations (RotD100/RotD50)
        """
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        rd100_resid_output = os.path.join(a_dstdir, "%s-%d-resid-rd100.txt" %
                                          (self.comp_label, sim_id))
        for comp in ['rotd50', 'rotd100', 'ratio']:
            # Build paths and check lengths
            fileroot = os.path.join(a_dstdir, "%s-%d_r0-%d-rd100-%s" %
                                    (self.comp_label, sim_id,
                                     self.max_cutoff, comp))
            bband_utils.check_path_lengths([rd100_resid_output, fileroot],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR,
                                         "resid2uncer_varN")) +
                   "residfile=%s fileroot=%s " %
                   (rd100_resid_output, fileroot) +
                   "comp=%s nstat=%d nper=63 " %
                   (comp, len(site_list)) +
                   "min_cdst=%d max_cdst=%d >> %s 2>&1" %
                   (0, self.max_cutoff, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Generate bias plot
        plot_mode = 'rd100'
        fileroot = ("%s-%d_r0-%d-rd100" %
                    (self.comp_label, sim_id,
                     self.max_cutoff))
        plottitle = ("GOF Comparison between %s and simulation %d" %
                     (self.comp_label, sim_id))
        plotter = PlotGoF()
        plotter.plot(plottitle, fileroot, a_dstdir, a_dstdir,
                     cutoff=self.max_cutoff, mode=plot_mode, colorset='single')

    def run(self):
        if self.obs_format is None:
            # simulation mode
            self.run_simulation()
        else:
            # validation mode
            self.run_validation()

    def run_simulation(self):
        """
        Generate RotD50/100 values for the calculated timeseries
        """
        print("RotDXX".center(80, '-'))
        #
        # convert input bbp acc files to peer format acc files
        #

        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.rotd100_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_outdir], print_cmd=False)

        self.calculate_simulated(a_statfile, a_tmpdir, a_outdir, a_outdir)

        # All done!
        print("RotDXX Completed".center(80, '-'))

    def run_validation(self):
        """
        Do all steps needed for creating the ratio of maximum to median
        response across orientations (RotD100/RotD50)
        """
        print("RotDXX".center(80, '-'))

        # Initialize
        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.rotd100_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir, "validations",
                                           "baker_rd100")

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_indir, a_tmpdir, a_tmpdir_seis,
                            a_outdir, a_validation_outdir],
                           print_cmd=False)

        # Source file, parse it!
        a_srcfile = os.path.join(a_indir, self.r_srcfile)
        self.src_keys = bband_utils.parse_src_file(a_srcfile)

        # Calculate RotD100/RotD50 for simulated seismograms
        self.calculate_simulated(a_statfile, a_tmpdir,
                                 a_outdir, a_validation_outdir)

        # Calculate RotD100/RotD50 for observation seismograms
        self.calculate_observations(a_indir, a_statfile,
                                    a_tmpdir_seis, a_validation_outdir)

        # Calculate ratios for simulated and observation data
        self.calculate_ratios(a_statfile, a_validation_outdir)

        # Generate comparison data table
        self.calculate_residuals(a_statfile, a_validation_outdir)

        # Generate bias plot showing the comparison between
        # simulations and observations
        self.generate_plot(a_statfile, a_validation_outdir)

        # All done!
        print("RotDXX Completed".center(80, '-'))

if __name__ == '__main__':
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    if len(sys.argv) == 3:
        # simulation mode
        ME = RotD100(sys.argv[1], sim_id=int(sys.argv[2]))
    else:
        # validation mode
        ME = RotD100(sys.argv[2], sys.argv[1],
                     sys.argv[3], sys.argv[4],
                     sys.argv[5], sys.argv[6],
                     sys.argv[7], sim_id=int(sys.argv[8]))
    ME.run()
