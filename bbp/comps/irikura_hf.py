#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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
import velocity_models
from station_list import StationList
from install_cfg import InstallCfg
from irikura_hf_cfg import IrikuraHFCfg, calculate_rvfac

class IrikuraHF(object):
    """
    This class contains the glue code needed to interface the
    Irikura Recipe Method 2 High Frequency and the Broadband Platform.
    """
    def __init__(self, i_r_srcfile, i_r_srffile, i_r_velmodel,
                 i_r_stations, vmodel_name, sim_id=0,
                 **kwargs):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_velmodel = i_r_velmodel
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.r_srcfiles = []
        self.stat_list = None
        self.install = None
        self.config = None
        self.log = None
        self.mean_rvfac = None
        self.range_rvfac = None

        # Get all src files that were passed to us
        if kwargs is not None and len(kwargs) > 0:
            for idx in range(len(kwargs)):
                self.r_srcfiles.append(kwargs['src%d' % (idx)])
        else:
            # Not a multisegment run, just use the single src file
            self.r_srcfiles.append(i_r_srcfile)

    def read_srf_header(self, a_srffile):
        """
        This function reads the SRF file header, returning a dictionary
        with the parameters found.
        """
        srf_dict = {}

        input_file = open(a_srffile, 'r')
        # Skip first two lines
        _ = input_file.readline()
        _ = input_file.readline()
        first_line_pieces = input_file.readline().strip().split()
        second_line_pieces = input_file.readline().strip().split()
        input_file.close()
        srf_dict["LON_TOP_CENTER"] = float(first_line_pieces[0])
        srf_dict["LAT_TOP_CENTER"] = float(first_line_pieces[1])
        srf_dict["WIDTH"] = float(first_line_pieces[5])
        srf_dict["STRIKE"] = float(second_line_pieces[0])
        srf_dict["DIP"] = float(second_line_pieces[1])
        srf_dict["DEPTH_TO_TOP"] = float(second_line_pieces[2])
        srf_dict["HYPO_ALONG_STK"] = float(second_line_pieces[3])
        srf_dict["HYPO_DOWN_DIP"] = float(second_line_pieces[4])

        return srf_dict

    def create_irikura_files(self, fault_param_dat_file, a_srffile):
        """
        This function creates the other files needed by the Irikura
        recipe after doing some math with the input parameters
        provided by the platform.
        """
        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        # Check for velocity model-specific parameters
        vmodel_params = vel_obj.get_codebase_params('gp')

        # Look for MEAN_RVFAC
        if 'MEAN_RVFAC' in vmodel_params:
            self.mean_rvfac = float(vmodel_params['MEAN_RVFAC'])
        else:
            self.mean_rvfac = self.config.VEL_RUP_FRAC
        # Look for RANGE_RVFAC
        if 'RANGE_RVFAC' in vmodel_params:
            self.range_rvfac = float(vmodel_params['RANGE_RVFAC'])
        else:
            self.range_rvfac = self.config.VEL_RUP_RANGE

        rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac, self.config.SEED)

        # Read SRF header
        srf_dict = self.read_srf_header(a_srffile)
        # Earth radius in km
        radius = 6371
        # Convert to radians
        tclat1 = math.pi / 180 * srf_dict["LAT_TOP_CENTER"]
        tclon1 = math.pi / 180 * srf_dict["LON_TOP_CENTER"]

        # Determine coodinate of center of fault plane
        azim = srf_dict["STRIKE"] + 90
        dist2 = srf_dict["WIDTH"] / 2 * math.cos(srf_dict["DIP"] *
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
        cent_dep = (srf_dict["DEPTH_TO_TOP"] * 1000 +
                    srf_dict["WIDTH"] / 2 * math.sin(srf_dict["DIP"] * math.pi /
                                                     180) * 1000) # In meters

        # Determine coordinate of hypocenter
        ddip = srf_dict["HYPO_DOWN_DIP"] * math.cos(srf_dict["DIP"] * math.pi /
                                                    180)
        # As HYPO_ALONG_STK --> 0.0, ddip/HYPO_ALONG_STK --> Inf,
        # atan(ddip/HYPO_ALONG_STK) --> pi/2
        if srf_dict["HYPO_ALONG_STK"] == 0.0:
            azim2 = srf_dict["STRIKE"] + math.pi / 2 * 180 / math.pi
        else:
            azim2 = (srf_dict["STRIKE"] +
                     math.atan(ddip / srf_dict["HYPO_ALONG_STK"]) *
                     180 / math.pi)
        dist3 = math.sqrt(srf_dict["HYPO_ALONG_STK"] ** 2 + ddip ** 2)

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
        hyp_dep = (srf_dict["DEPTH_TO_TOP"] * 1000 +
                   srf_dict["HYPO_DOWN_DIP"] * 1000 *
                   math.sin(srf_dict["DIP"] * math.pi / 180))

        # Determine density/Vs at "fault" and at "bedrock" (defined as
        # top of fault plane for now)

        # modified by A.I. 2019.7 =======
        #   depth of "ELEMENT FAULT" = cent_dep
        #   d_goto = Depth of "SEISMIC-BEDROCK"

        vs_m = self.config.vmodel["vs_m"]
        rho = self.config.vmodel["rho"]

        # ELEMENT-FAULT
        index = bisect.bisect_left(self.config.vmodel["depth0"], cent_dep)
        density = rho[index] # g/cm3

        # density_bed = density
        fault_vs = float(vs_m[index]) / 1000
        # fault_vs_bed = fault_vs

        # Depth of SEISMIC-BEDROCK
        if srf_dict["DEPTH_TO_TOP"] < 1.0:
            d_goto = 1000
        else:
            d_goto = srf_dict["DEPTH_TO_TOP"] * 1000

        index_bed = bisect.bisect_left(self.config.vmodel["depth0"], d_goto)
        density_bed = rho[index_bed] # g/cm3
        fault_vs_bed = float(vs_m[index_bed]) / 1000
        # ==== A.I. 2019.7 END

        dlen = self.config.DXX / 1000.0
        dwid = self.config.DYY / 1000.0

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
        desvel = vs_pq * rvfac

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
        out_file.write("100.0 \t\t# SAMPLING-FREQUENCY(Hz)\n")
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
        out_file.write("%.3f \t\t# L\n" % (dlen))
        out_file.write("%.3f \t\t# W\n" % (dwid))
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
        for idx, _ in enumerate(self.stat_list.get_station_list(), 1):
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
        for idx, _ in enumerate(self.stat_list.get_station_list(), 1):
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
        for idx, station in enumerate(self.stat_list.get_station_list(), 1):
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

        for idx, station in enumerate(self.stat_list.get_station_list()):
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
        a_param_outdir = os.path.join(a_outdir, "param_files")

        #
        # Make sure the output and two tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, irikura_dir, a_outdir, a_param_outdir,
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
        a_sdropout = os.path.join(a_tmpdir,
                                  self.config.sdropout)
        a_segments_midpoint = os.path.join(a_tmpdir,
                                           self.config.segments_midpoint)
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
        elem_param_rndt_dat_file = os.path.join(irikura_dir,
                                                "%s_elem_param_rndt.dat" %
                                                (sim_id))

        # Copy needed files to Irikura Directory
        #  - velocity model
        #  - segments_midpoint
        #  - sdropout
        #  - srf file(s)
        shutil.copy2(a_velmodel,
                     os.path.join(irikura_dir, self.r_velmodel))
        shutil.copy2(a_sdropout,
                     os.path.join(irikura_dir, self.config.sdropout))
        shutil.copy2(a_segments_midpoint,
                     os.path.join(irikura_dir, self.config.segments_midpoint))
        shutil.copy2(a_srffile,
                     os.path.join(irikura_dir, self.r_srffile))
        if len(self.r_srcfiles) == 1:
            # No need to copy anything else, set up filenames for all SRF files
            r_single_seg_srf = self.r_srffile
        else:
            # Copy all SRF files
            r_single_seg_srf = "single_seg.%s" % (self.r_srffile)
            a_single_seg_srf = os.path.join(a_tmpdir, r_single_seg_srf)
            shutil.copy2(a_single_seg_srf,
                         os.path.join(irikura_dir, r_single_seg_srf))
            for segnum in range(len(self.r_srcfiles)):
                r_seg_srf = "seg%d.%s" % ((segnum + 1), self.r_srffile)
                a_seg_srf = os.path.join(a_tmpdir, r_seg_srf)
                shutil.copy2(a_seg_srf,
                             os.path.join(irikura_dir, r_seg_srf))

        # Create Irikura velocity model file
        self.create_velocity_file(vel_file, vel_file_p)
        shutil.copy2(vel_file,
                     os.path.join(a_param_outdir, "%d_soil.dat" % (sim_id)))
        shutil.copy2(vel_file_p,
                     os.path.join(a_param_outdir, "%d_soil_p.dat" % (sim_id)))

        # Create Irikura station list
        self.create_station_list(station_file)
        shutil.copy2(station_file,
                     os.path.join(a_param_outdir, "%d_station.dat" % (sim_id)))

        # Create other input files
        self.create_irikura_files(fault_param_dat_file,
                                  os.path.join(irikura_dir, r_single_seg_srf))
        shutil.copy2(fault_param_dat_file,
                     os.path.join(a_param_outdir, "%d_fault_param.dat" % (sim_id)))

        # Create phase files
        r_phase_file = "phase.dat"
        r_phase2_file = "phase2.dat"
        a_phase_file = os.path.join(irikura_dir, r_phase_file)
        a_phase2_file = os.path.join(irikura_dir, r_phase2_file)
        self.create_phase_files(a_phase_file, a_phase2_file)

        # Save a copy in param_dir
        shutil.copy2(a_phase_file, os.path.join(a_param_outdir, r_phase_file))
        shutil.copy2(a_phase2_file, os.path.join(a_param_outdir, r_phase2_file))

        # Irikura binaries
        srf2grns_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                    "srf2grns")
        dtrandom_bin = os.path.join(install.A_IRIKURA_BIN_DIR,
                                    "dtrandom")
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

        # Mesh size parameters
        dxx = self.config.DXX
        dyy = self.config.DYY

        # Open combined elem_param.dat file
        elem_param_file = open(elem_param_dat_file, 'w')
        istart = 0

        for segnum in range(len(self.r_srcfiles)):
            if len(self.r_srcfiles) == 1:
                # Single segment simulation
                r_seg_srf = self.r_srffile
            else:
                # Multi segment simulation
                r_seg_srf = "seg%d.%s" % ((segnum + 1), self.r_srffile)

            current_elem_param_dat_file = os.path.join(irikura_dir,
                                                        "%s_elem_param_seg%d.dat" %
                                                        (sim_id, segnum + 1))
            srf2grns_input_file = os.path.join(irikura_dir,
                                               "input%d.txt" % (segnum + 1))

            # Write config file for srf2grns
            config_file = open(srf2grns_input_file, 'w')
            config_file.write("%d\n" % (segnum + 1))
            config_file.write("%s\n" % (r_single_seg_srf))
            config_file.write("%s\n" % (r_seg_srf))
            config_file.write("%s\n" % (self.config.segments_midpoint))
            config_file.write("%s\n" % (self.config.sdropout))
            config_file.write("%s\n" % (os.path.basename(current_elem_param_dat_file)))
            config_file.write("%f %f\n" % (dxx, dyy))
            config_file.close()

            # Run the srf2grns code
            cmd = ("%s < %s >> %s 2>&1" %
                   (srf2grns_bin,
                    os.path.basename(srf2grns_input_file),
                    self.log))
            bband_utils.runprog(cmd, abort_on_error=True)

            # Save copy
            shutil.copy2(current_elem_param_dat_file,
                         os.path.join(a_param_outdir,
                                      os.path.basename(current_elem_param_dat_file)))

            current_elem_file = open(current_elem_param_dat_file, 'r')
            for line in current_elem_file:
                line = line.strip()
                if not line:
                    continue
                pieces = line.split()
                pieces[0] = int(float(pieces[0]))
                pieces[0] = pieces[0] + istart
                elem_param_file.write("%d %s %s %s %s %s %s\n" % (pieces[0],
                                                                  pieces[1],
                                                                  pieces[2],
                                                                  pieces[3],
                                                                  pieces[4],
                                                                  pieces[5],
                                                                  pieces[6]))
            current_elem_file.close()
            istart = pieces[0]

        elem_param_file.close()

        # Save copy
        shutil.copy2(elem_param_dat_file,
                     os.path.join(a_param_outdir,
                                  os.path.basename(elem_param_dat_file)))

        # Now add rupture perturbations by using the dtrandom code
        dtrandom_input_file = os.path.join(irikura_dir, "dtrandom_input.txt")
        dtrdm = open(dtrandom_input_file, 'w')
        dtrdm.write("%d\n" % (self.config.SEED))
        dtrdm.write("%d\n" % (istart))
        elem_param_file = open(elem_param_dat_file, 'r')
        for line in elem_param_file:
            line = line.strip()
            if not line:
                continue
            dtrdm.write("%s\n" % (line))
        elem_param_file.close()
        dtrdm.close()

        # Save copy
        shutil.copy2(dtrandom_input_file, os.path.join(a_param_outdir,
                                                       "dtrandom_input.txt"))

        # Run the dtrandom code
        cmd = ("%s < %s > %s" %
               (dtrandom_bin,
                dtrandom_input_file,
                elem_param_rndt_dat_file))
        bband_utils.runprog(cmd, abort_on_error=True)

        # Save copy
        shutil.copy2(elem_param_rndt_dat_file,
                     os.path.join(a_param_outdir,
                                  os.path.basename(elem_param_rndt_dat_file)))

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
                   (greenscale_bin, station_file, elem_param_rndt_dat_file,
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
