#!/usr/bin/env python
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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import math
import bisect
import random
import shutil

# Import Broadband modules
import bband_utils
import validation_cfg
from station_list import StationList
from install_cfg import InstallCfg
from irikura_hf_cfg import IrikuraHFCfg

class IrikuraHF(object):
    """
    This class contains the glue code needed to interface the
    Irikura Recipe Method 2 High Frequency and the Broadband Platform.
    """
    def __init__(self, i_r_srcfile, i_r_srffile, i_r_velmodel,
                 i_r_stations, vmodel_name,
                 val_name=None, sim_id=0):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_velmodel = i_r_velmodel
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.val_name = val_name
        self.val_obj = None
        self.stat_list = None
        self.install = None
        self.config = None
        self.log = None

    def create_irikura_files(self, fault_param_dat_file):
        """
        This function creates the other files needed by the Irikura
        recipe after doing some math with the input parameters
        provided by the platform.
        """
        # Earth radius in km
        radius = 6371
        # Convert to radians
        tclat1 = math.pi / 180 * self.config.LAT_TOP_CENTER
        tclon1 = math.pi / 180 * self.config.LON_TOP_CENTER

        # Determine coodinate of center of fault plane
        azim = self.config.STRIKE + 90
        dist2 = self.config.WIDTH / 2 * math.cos(self.config.DIP *
                                                 math.pi / 180)

        cent_lat = 180 / math.pi * (math.asin(math.sin(tclat1) *
                                              math.cos(dist2/radius) +
                                              math.cos(tclat1) *
                                              math.sin(dist2/radius) *
                                              math.cos(azim * math.pi / 180)))
        cent_lon = 180 / math.pi * (tclon1 +
                                    math.atan2(math.sin(azim * math.pi / 180) *
                                               math.sin(dist2/radius) *
                                               math.cos(tclat1),
                                               math.cos(dist2/radius) -
                                               math.sin(tclat1) *
                                               math.sin(cent_lat *
                                                        math.pi / 180)))
        cent_dep = (self.config.DEPTH_TO_TOP * 1000 +
                    self.config.WIDTH / 2 * math.sin(self.config.DIP * math.pi /
                                                     180) * 1000) # In meters

        # Determine coordinate of hypocenter
        ddip = self.config.HYPO_DOWN_DIP * math.cos(self.config.DIP * math.pi /
                                                    180)
        # As HYPO_ALONG_STK --> 0.0, ddip/HYPO_ALONG_STK --> Inf,
        # atan(ddip/HYPO_ALONG_STK) --> pi/2
        if self.config.HYPO_ALONG_STK == 0.0:
            azim2 = self.config.STRIKE + math.pi / 2 * 180 / math.pi
        else:
            azim2 = (self.config.STRIKE +
                     math.atan(ddip / self.config.HYPO_ALONG_STK) *
                     180 / math.pi)
        dist3 = math.sqrt(self.config.HYPO_ALONG_STK ** 2 + ddip ** 2)

        hyp_lat = 180 / math.pi * (math.asin(math.sin(tclat1) *
                                             math.cos(dist3/radius) +
                                             math.cos(tclat1) *
                                             math.sin(dist3/radius) *
                                             math.cos(azim2 * math.pi / 180)))
        hyp_lon = 180 / math.pi * (tclon1 +
                                   math.atan2(math.sin(azim2 * math.pi / 180) *
                                              math.sin(dist3/radius) *
                                              math.cos(tclat1),
                                              math.cos(dist3/radius) -
                                              math.sin(tclat1) *
                                              math.sin(hyp_lat *
                                                       math.pi / 180)))
        hyp_dep = (self.config.DEPTH_TO_TOP * 1000 +
                   self.config.HYPO_DOWN_DIP * 1000 *
                   math.sin(self.config.DIP * math.pi / 180))

        # Determine density/Vs at "fault" and at "bedrock" (defined as
        # top of fault plane for now)
        d_goto = cent_dep

        # As per Pylint, this variable is not used
        # depth = self.config.vmodel["depth"]
        vs_m = self.config.vmodel["vs_m"]
        rho = self.config.vmodel["rho"]

        index = bisect.bisect_left(self.config.vmodel["depth0"], d_goto)
        density = rho[index] # g/cm3
        density_bed = density
        fault_vs = float(vs_m[index]) / 1000
        fault_vs_bed = fault_vs

        dwid = 2
        dlen = 2

        # AI: stress drop on element = 10MPa 2018.7.25
        sigma_s = 1.0e7  # Pa

        # AI: mo_elem modified 2018.7.25
        r_elem = math.sqrt(dwid * dlen * 1.0e6 / math.pi)
        mo_elem = sigma_s * (16. / 7.) * (r_elem **3)  # Nm

        # Units g/cm3
        ro_pq = density
        # Units km/s
        vs_pq = fault_vs

        myu = ro_pq * vs_pq * vs_pq * (1000 ** 3) # Units Pa
        ds_m = mo_elem / myu / dwid / dlen / 1000 / 1000 # Units m

        # Units km/s
        desvel = vs_pq * 0.72

        # AI: radiation pattern 0.63/sqrt(2) and m_elem 2018.7.25
        #radpatt = 0.4384
        radpatt = 0.445
        m_elem = 4

        # Output the DAT file
        out_file = open(fault_param_dat_file, 'w')
        out_file.write("# FAULT-PARAM.DAT\n")
        out_file.write("%.6f \t# FAULT-PLANE-CENTER-LATITUDE\n" % (cent_lat))
        out_file.write("%.6f \t# FAULT-PLANE-CENTER-LONGITUDE\n" % (cent_lon))
        out_file.write("%.4f \t# FAULT-PLANE-CENTER-DEPTH(GL-m)\n" % (cent_dep))
        out_file.write("%.8f \t\t# HYPOCENTER-LATITUDE\n" % (hyp_lat))
        out_file.write("%.8f \t\t# HYPOCENTER-LONGITUDE\n" % (hyp_lon))
        out_file.write("%.4f \t\t# HYPOCENTER-DEPTH(GL-m)\n" % (hyp_dep))
        out_file.write("120.0 \t\t# SAMPLING-FREQUENCY(Hz)\n")
        out_file.write("120.0 \t\t# SAMPLING-TIME(sec)\n")
        out_file.write("%.2f \t\t# ELEMENT-DENSITY(g/cm3)\n" % (density))
        out_file.write("%.2f \t\t# SEISMIC-BEDROCK-DENSITY(g/cm3)\n" %
                       (density_bed))
        out_file.write("%.2f \t\t# ELEMENT-SWAVE-VELOCITY(km/sec)\n" %
                       (fault_vs))
        out_file.write("%.2f \t\t# SEISMIC-BEDROCK-SWAVE-VELOCITY(km/sec)\n" %
                       (fault_vs_bed))
        mo_string = "%8.2E \t# ELEMENT-MOMENT(M0s)\n" % (mo_elem)
        mo_string = mo_string.replace('E+0', 'E+')
        out_file.write("%s" % (mo_string))
        stress_string = "%.2E \t# ELEMENT-STRESS-DROP(PA)\n" % (sigma_s)
        stress_string = stress_string.replace('E+0', 'E+')
        out_file.write("%s" % (stress_string))
        out_file.write("%.4f \t\t# RADIATION-PATTERN\n" % (radpatt))
        out_file.write("6.0 \t\t# FMAX\n")
        out_file.write("%.3f \t\t# M\n" % (m_elem))
        out_file.write("100.0 \t\t# Q-VALUE(A1) "
                       "Q = A1 * f(Hz)^A2(A3 < f(Hz))\n")
        out_file.write("0.69 \t\t# Q-VALUE(A2) Q = A4( A3 > f(Hz) )\n")
        out_file.write("1.0 \t\t# Q-VALUE(A3)\n")
        out_file.write("100.0 \t\t# Q-VALUE(A4)\n")
        out_file.write("# green scale param\n")
        out_file.write("2.0 \t\t# L\n")
        out_file.write("2.0 \t\t# W\n")
        out_file.write("%f \t\t# Ds\n" % (ds_m))
        myu_string = "%.2E \t\t# myu\n" % (myu)
        myu_string = myu_string.replace('E+0', 'E+')
        out_file.write("%s" % (myu_string))
        out_file.write("%.3f \t\t# Vr\n" % (desvel))
        out_file.write("0.0 \t\t# offset\n")
        out_file.close()

    def create_velocity_file(self, vel_file, vel_file_p):
        """
        This function creates the Irikura velocity model file
        """
        # Get the parameters we need
        thick = self.config.vmodel["h"]
        vs_km = self.config.vmodel["vs"]
        vp_km = self.config.vmodel["vp"]
        rho = self.config.vmodel["rho"]
        qs = self.config.vmodel["qs"]
        # Convert to meters
        thick = [item * 1000 for item in thick]
        vs_m = [item * 1000 for item in vs_km]
        vp_m = [item * 1000 for item in vp_km]

        # Convert depths to absolute depths
        depth = [0] * len(thick)
        depth[0] = thick[0]
        for idx in range(1, len(thick)):
            depth[idx] = depth[idx - 1] + thick[idx]
        depth[len(thick) - 1] = 9999999999
        # Make a copy of the original array
        self.config.vmodel["depth0"] = depth[:]

        if self.config.DEPTH_TO_TOP < 1.0:
            idx = bisect.bisect_left(depth, 1000)
            if len(depth) == idx:
                raise bband_utils.ParameterError("Velocity model above ztor!")
            depth = depth[:idx+1]
        else:
            # Now, only pick the ones up to fault depth
            idx = bisect.bisect_left(depth, self.config.DEPTH_TO_TOP * 1000)
            if len(depth) == idx:
                raise bband_utils.ParameterError("Velocity model above ztor!")
            if not depth[idx] == self.config.DEPTH_TO_TOP * 1000:
                depth[idx] = self.config.DEPTH_TO_TOP * 1000
            # Select values up to the top of the fault
            depth = depth[:idx+1]

        # Write vel_file
        out_file = open(vel_file, 'w')
        out_file.write("# SOIL-LAYER.DAT\n")
        out_file.write("# SOIL-PARAM(Vs,Ro,Qs)\n")
        # Write Vs, Rho, Qs
        for idx in range(0, len(depth) + 1):
            out_file.write("%d %d %4.2f %d\n" %
                           (idx + 1, vs_m[idx], rho[idx], round(qs[idx])))
        out_file.write("# LAYER-DEPTH(GL-m)\n")
        for idx, _ in enumerate(self.stat_list.getStationList(), 1):
            out_file.write("%d %d %d" % (idx, 0, 0))
            for item in depth:
                out_file.write(" %.1f" % (item))
            out_file.write("\n")
        out_file.close()

        # Write vel_file_p
        out_file = open(vel_file_p, 'w')
        out_file.write("# SOIL-LAYER.DAT\n")
        out_file.write("# SOIL-PARAM(Vp,Ro,Qs)\n")
        # Write Vs, Rho, Qs
        for idx in range(0, len(depth) + 1):
            out_file.write("%d %d %4.2f %d\n" %
                           (idx + 1, vp_m[idx], rho[idx], round(qs[idx])))
        out_file.write("# LAYER-DEPTH(GL-m)\n")
        for idx, _ in enumerate(self.stat_list.getStationList(), 1):
            out_file.write("%d %d %d" % (idx, 0, 0))
            for item in depth:
                out_file.write(" %.1f" % (item))
            out_file.write("\n")
        out_file.close()

        # Store velocity model results
        self.config.vmodel["depth"] = depth
        self.config.vmodel["vs_m"] = vs_m
        self.config.vmodel["vp_m"] = vp_m

    def create_station_list(self, station_file):
        """
        This function creates the Irikura recipe station list
        """
        out_file = open(station_file, 'w')
        for idx, station in enumerate(self.stat_list.getStationList(), 1):
            out_file.write("%d\t%9.5f\t%9.5f\t0.0  \tSAS%s.dat\n" %
                           (idx, station.lat, station.lon, station.scode))
        out_file.close()

    def create_phase_files(self, phase_file, phase2_file):
        """
        This function creates the two phase.dat files using
        random numbers between [-pi:pi]
        """
        output_filenames = [phase_file, phase2_file]
        values_needed = [16386, 16386]

        for output_filename, num_values in zip(output_filenames,
                                               values_needed):
            output_file = open(output_filename, 'w')
            for _ in range(0, num_values):
                output_file.write("%2.7f\n" %
                                  (random.random() * 2 * math.pi - math.pi))
            output_file.close()

    def convert_to_bbp(self, in_hor_file, in_ver_file, out_acc_file):
        """
        This function converts the in_file Irikura seismogram to an
        acceleration BBP file
        """
        header_lines = 0
        freq = None
        fac1_hor = None
        fac2_hor = None
        fac1_ver = None
        fac2_ver = None

        # Read horizontal acc file
        irikura_file = open(in_hor_file, 'r')
        for line in irikura_file:
            line = line.strip()
            header_lines = header_lines + 1
            if line.startswith("Sampling Freq"):
                pieces = line.split()
                token = pieces[2]
                freq = float(token[0:token.find("Hz")])
                continue
            if line.startswith("Scale Factor"):
                pieces = line.split()
                token = pieces[2]
                fac2_hor = float(token.split("/")[1])
                fac1_hor = float(token.split("/")[0][0:token.find("(gal)")])
                continue
            if line.startswith("Memo"):
                break
        irikura_file.close()

        # Read vertical acc file
        irikura_file = open(in_ver_file, 'r')
        for line in irikura_file:
            line = line.strip()
            if line.startswith("Scale Factor"):
                pieces = line.split()
                token = pieces[2]
                fac2_ver = float(token.split("/")[1])
                fac1_ver = float(token.split("/")[0][0:token.find("(gal)")])
                continue
            if line.startswith("Memo"):
                break
        irikura_file.close()

        if (freq is None
            or fac1_hor is None
            or fac2_hor is None
            or fac1_ver is None
            or fac2_ver is None):
            # Not able to parse it properly, exit!
            raise bband_utils.ProcessingError("Could not parse files %s and %s" %
                                              (in_hor_file, in_ver_file))

        time_step = 1.0 / freq
        time_curr = 0.0
        skip_lines = 0
        irikura_hor_file = open(in_hor_file, 'r')
        irikura_ver_file = open(in_ver_file, 'r')
        bbp_file = open(out_acc_file, 'w')
        # Add header to BBP file
        bbp_file.write("#    time(sec)      N-S(cm/s/s)      "
                       "E-W(cm/s/s)      U-D(cm/s/s)\n")
        for line1, line2 in zip(irikura_hor_file, irikura_ver_file):
            # Skip header lines
            if skip_lines < header_lines:
                skip_lines = skip_lines + 1
                continue
            line1 = line1.strip()
            line2 = line2.strip()
            pieces1 = line1.split()
            pieces1 = [float(val) for val in pieces1]
            pieces2 = line2.split()
            pieces2 = [float(val) for val in pieces2]
            # Convert each value to gal (cm/s/s)
            pieces1 = [val * fac1_hor / fac2_hor for val in pieces1]
            pieces2 = [val * fac1_ver / fac2_ver for val in pieces2]
            # Write values to BBP file, repeating as needed for the
            # 3-component file
            for piece_h, piece_v in zip(pieces1, pieces2):
                bbp_file.write("%15.6e%15.6e%15.6e%15.6e\n" %
                               (time_curr, piece_h, piece_h, piece_v))
                # Don't forget to increment time
                time_curr = time_curr + time_step
        bbp_file.close()
        irikura_hor_file.close()
        irikura_ver_file.close()

    def create_vel_bbp(self, stat):
        """
        This function derives a velocity bbp file from an acceleration
        file
        """
        install = self.install
        sim_id = self.sim_id
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        in_acc_file = os.path.join(a_tmpdir, "%d.%s-hf.acc.bbp" %
                                   (sim_id, stat))
        out_vel_file = os.path.join(a_tmpdir, "%d.%s-hf.bbp" %
                                    (sim_id, stat))

        # Since we have acceleration files, we need to integrate to
        # get to velocity

        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.acc.000" % (sim_id, stat))
        ewfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.acc.090" % (sim_id, stat))
        udfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.acc.ver" % (sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
               " wcc2bbp=0 nsfile=%s ewfile=%s udfile=%s < %s >> %s 2>&1" %
               (nsfile, ewfile, udfile, in_acc_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True)
        # Now we need to integrate to get to velocity
        for comp in ['000', '090', 'ver']:
            file_in = os.path.join(a_tmpdir,
                                   "%d.%s-hf.acc.%s" % (sim_id, stat, comp))
            file_out = os.path.join(a_tmpdir,
                                    "%s.%s-hf.vel.%s" % (sim_id, stat, comp))
            bband_utils.check_path_lengths([file_in, file_out],
                                           bband_utils.GP_MAX_FILENAME)
            cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                   " integ=1 filein=%s fileout=%s" % (file_in, file_out))
            bband_utils.runprog(cmd, abort_on_error=True)
        # Now we put together the file again as a velocity bbp file
        nsfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.vel.000" % (sim_id, stat))
        ewfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.vel.090" % (sim_id, stat))
        udfile = os.path.join(a_tmpdir,
                              "%d.%s-hf.vel.ver" % (sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
               " wcc2bbp=1 nsfile=%s ewfile=%s udfile=%s > %s 2>> %s" %
               (nsfile, ewfile, udfile, out_vel_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True)

    def process_seismograms(self, irikura_dir):
        """
        This function reads the seismograms generated by the Irikura
        method and generates the velocity and acceleration
        seismograms in the BBP format needed by the platform
        """
        acc_hor_dir = os.path.join(irikura_dir, "HOR", "acc")
        acc_ver_dir = os.path.join(irikura_dir, "VER", "acc")
        tmp_dir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))

        for idx, station in enumerate(self.stat_list.getStationList()):
            in_hor_file = os.path.join(acc_hor_dir, "ac%06d.dat" % (idx + 1))
            in_ver_file = os.path.join(acc_ver_dir, "ac%06d.dat" % (idx + 1))
            out_acc_file = os.path.join(tmp_dir, "%d.%s-hf.acc.bbp" %
                                        (self.sim_id, station.scode))
            self.convert_to_bbp(in_hor_file, in_ver_file, out_acc_file)
            self.create_vel_bbp(station.scode)

    def run(self):
        """
        This function prepares the parameter file for the Irikura
        Recipe, invokes it, and formats its output to be compatible
        with the Broadband Platform
        """
        print("IrikuraHF".center(80, '-'))

        self.install = InstallCfg.getInstance()
        install = self.install
        sim_id = self.sim_id

        # Find validation object if this is a validation run
        if self.val_name is not None:
            self.val_obj = validation_cfg.VE_EVENTS.get_event_by_name(self.val_name)

        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.irikura_hf_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        irikura_dir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                   "irikura_hf_%s" % (sta_base))
        irikura_hor_dir = os.path.join(irikura_dir, "HOR")
        irikura_ver_dir = os.path.join(irikura_dir, "VER")
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        #
        # Make sure the output and two tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, irikura_dir, a_outdir,
                            irikura_hor_dir, irikura_ver_dir])

        a_velmodel = os.path.join(a_indir, self.r_velmodel)
        a_stations = os.path.join(a_indir, self.r_stations)
        a_srffile = os.path.join(a_indir, self.r_srffile)
        self.stat_list = StationList(a_stations)

        # Read input files, calculate CSM parameters
        self.config = IrikuraHFCfg(os.path.join(install.A_IN_DATA_DIR,
                                                str(sim_id),
                                                self.r_srcfile),
                                   a_velmodel)

        # Initilize random number generator with seed from src file
        random.seed(self.config.SEED)

        # Create filenames for all intermediate files
        vel_file = os.path.join(irikura_dir,
                                "%d_soil.dat" %
                                (sim_id))
        vel_file_p = os.path.join(irikura_dir,
                                  "%d_soil_p.dat" %
                                  (sim_id))
        station_file = os.path.join(irikura_dir,
                                    "%d_station.dat" %
                                    (sim_id))
        fault_param_dat_file = os.path.join(irikura_dir,
                                            "%d_fault_param.dat" %
                                            (sim_id))
        elem_param_dat_file = os.path.join(irikura_dir,
                                           "%s_elem_param.dat" %
                                           (sim_id))
        srf2grns_input_file = os.path.join(irikura_dir,
                                           "input.txt")
        phase_file = os.path.join(irikura_dir, "phase.dat")
        phase2_file = os.path.join(irikura_dir, "phase2.dat")

        # Create Irikura velocity model file
        self.create_velocity_file(vel_file, vel_file_p)

        # Create Irikura station list
        self.create_station_list(station_file)

        # Create other input files
        self.create_irikura_files(fault_param_dat_file)

        # Create phase files
        self.create_phase_files(phase_file, phase2_file)

        # Copy velocity model and srf file to Irikura dir
        shutil.copy2(a_velmodel,
                     os.path.join(irikura_dir, self.r_velmodel))
        shutil.copy2(a_srffile,
                     os.path.join(irikura_dir, self.r_srffile))

        # Irikura binaries
        srf2grns_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                    "srf2grns")
        statgreen_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                     "statgreen")
        statgreen2_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                      "statgreen2")
        greenscale_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                      "greenscale")

        # Run in tmpdir subdir to isolate files
        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(irikura_dir)

        # Run the srf2grns code
        config_file = open(srf2grns_input_file, 'w')
        config_file.write("%s\n" % (self.r_srffile))
        config_file.write("%s\n" % (os.path.basename(elem_param_dat_file)))
        config_file.close()
        cmd = ("%s < %s >> %s 2>&1" %
               (srf2grns_bin,
                os.path.basename(srf2grns_input_file),
                self.log))

        bband_utils.runprog(cmd, abort_on_error=True)

        # Run the statgreen program
        os.chdir(irikura_hor_dir)
        cmd = ("%s %s %s %s ../phase.dat >> %s 2>&1" % (statgreen_bin,
                                                        station_file,
                                                        fault_param_dat_file,
                                                        vel_file,
                                                        self.log))
        bband_utils.runprog(cmd, abort_on_error=True)

        # Create directories for the AS and SAS files, and move data
        bband_utils.mkdirs([os.path.join(irikura_hor_dir, "AS"),
                            os.path.join(irikura_hor_dir, "SAS")])
        for fname in glob.iglob(os.path.join(irikura_hor_dir,
                                             "AS*.dat")):
            shutil.move(fname, os.path.join(irikura_hor_dir, "AS"))
        for fname in glob.iglob(os.path.join(irikura_hor_dir,
                                             "SAS*.dat")):
            shutil.move(fname, os.path.join(irikura_hor_dir, "SAS"))

        # Run statgreen2 (for vertical component)
        os.chdir(irikura_ver_dir)
        cmd = ("%s %s %s %s ../phase2.dat >> %s 2>&1" % (statgreen2_bin,
                                                         station_file,
                                                         fault_param_dat_file,
                                                         vel_file_p,
                                                         self.log))
        bband_utils.runprog(cmd, abort_on_error=True)

        # Create directories for the AS and SAS files, and move data
        bband_utils.mkdirs([os.path.join(irikura_ver_dir, "AS"),
                            os.path.join(irikura_ver_dir, "SAS")])
        for fname in glob.iglob(os.path.join(irikura_ver_dir,
                                             "AS*.dat")):
            shutil.move(fname, os.path.join(irikura_ver_dir, "AS"))
        for fname in glob.iglob(os.path.join(irikura_ver_dir,
                                             "SAS*.dat")):
            shutil.move(fname, os.path.join(irikura_ver_dir, "SAS"))

        # Now, run the greenscale code for the
        # horizontal and vertical components
        for working_dir in [irikura_hor_dir, irikura_ver_dir]:
            os.chdir(working_dir)

            cmd = ("%s SAS %s %s %s LogGS.dat HP 1.0 0 >> %s 2>&1" %
                   (greenscale_bin, station_file, elem_param_dat_file,
                    fault_param_dat_file, self.log))
            bband_utils.runprog(cmd, abort_on_error=True)

            # Create directories for the output files, move data
            bband_utils.mkdirs([os.path.join(working_dir, "acc"),
                                os.path.join(working_dir, "vel"),
                                os.path.join(working_dir, "velf")])
            for fname in glob.iglob(os.path.join(working_dir,
                                                 "ac0*.dat")):
                shutil.move(fname, os.path.join(working_dir, "acc"))

            for fname in glob.iglob(os.path.join(working_dir,
                                                 "ve0*HP*.dat")):
                shutil.move(fname, os.path.join(working_dir, "velf"))
            for fname in glob.iglob(os.path.join(working_dir,
                                                 "ve0*.dat")):
                shutil.move(fname, os.path.join(working_dir, "vel"))

        # Restore working directory
        os.chdir(old_cwd)

        # Need to copy and re-format output seismograms
        self.process_seismograms(irikura_dir)

        print("IrikuraHF Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (sys.argv[0]))
    ME = IrikuraHF(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                   sys.argv[5], sys.argv[6], sim_id=int(sys.argv[7]))
    ME.run()
