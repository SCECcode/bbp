#!/bin/env python3
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

ExSim Broadband Platform module
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import velocity_models
from station_list import StationList
from install_cfg import InstallCfg
from exsim_cfg import ExSimCfg

class ExSim(object):
    """
    This class contains the glue code needed to interface ExSim and
    the Broadband Platform.
    """
    def __init__(self, i_r_srcfile, i_r_param_template,
                 i_r_stations, vmodel_name, sim_id=0):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_param_template = i_r_param_template
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.stat_list = None
        self.num_stations = None
        self.install = None
        self.config = None
        self.kappa = None
        self.stress = None
        self.crustal_amp = None
        self.site_amp = None
        self.empirical_dir = None
        self.empirical_ranges = None
        self.slip_weights = None
        self.template_in = None
        self.param_out = None
        self.log = None

    def template_replace(self, key_list, value_list):
        """
        This function looks for the next line in the template input
        file containing all keys in key_list. It then replaces these
        keys with the values in value_list. Lines not matching keys
        are copied from the template to the output file. The
        operations stops as soon as a line matching all keys is found
        the replaced.
        """
        for line in self.template_in:
            copied = False
            for key in key_list:
                if key in line:
                    continue
                # Not the line we want, just copy it
                self.param_out.write(line)
                copied = True
                break
            if copied:
                # Copied already, look for the next line
                continue
            # Ok, this is the line we want, make substitutions
            for key, value in zip(key_list, value_list):
                line = line.replace(key, str(value), 1)
            # Write output line
            self.param_out.write(line)
            # All done!
            break

    def convert_exsim_to_rd50(self, in_exsim_file, out_rotd50):
        """
        This function converts an exsim output file to the rotd50
        format
        """
        # Open input and output files
        in_file = open(in_exsim_file, 'r')
        out_file = open(out_rotd50, 'w')

        # Skip input file header
        for line in in_file:
            if line.find("data format") >= 0:
                break

        # Write output header
        out_file.write("#  RotD50 RotD50 RotD50\n")
        out_file.write("#  %s\n" % (os.path.basename(in_exsim_file)))
        out_file.write("#  %s\n" % (os.path.basename(in_exsim_file)))
        out_file.write("#     63    0.0500\n")

        # Now convert the file to rotD50 format
        for line in in_file:
            line = line.strip()
            tokens = line.split()
            if len(tokens) != 6:
                # Invalid line
                continue
            period = tokens[1]
            value = tokens[3]
            try:
                period = float(period)
                value = float(value)
            except ValueError:
                # Not a valid line either
                continue
            # Convert value from cm/s/s to g
            value = value / bband_utils.G2CMSS
            # Write output
            out_file.write("%10.4f    %10.5e    %10.5e    %10.5e\n" %
                           (period, value, value, value))

        # Close files
        in_file.close()
        out_file.close()

    def convert_exsim_to_bbp(self, in_exsim_acc_file, stat,
                             tmpdir, outdir):
        """
        This function reads the in_exsim_acc_file, writes a bbp acc
        file and also generates the corresponding bbp vel file (both
        of those in the outdata directory)
        """
        out_bbp_acc_file = os.path.join(outdir, "%s.%s.acc.bbp" %
                                        (str(self.sim_id), stat))
        out_bbp_vel_file = os.path.join(outdir, "%s.%s.vel.bbp" %
                                        (str(self.sim_id), stat))
        # First write the BBP ACC file
        in_seis = open(in_exsim_acc_file, 'r')
        out_acc = open(out_bbp_acc_file, 'w')
        out_acc.write("# Sim stat=%s\n" % (stat))
        out_acc.write("#    time(sec)      N-S(cm/s/s)      "
                      "E-W(cm/s/s)      U-D(cm/s/s)\n")
        ready_to_copy = False
        for line in in_seis:
            if not ready_to_copy:
                line = line.strip()
                if line.find("time(s) acc(cm/s**2)") < 0:
                    # Not the line we are looking for...
                    continue
                ready_to_copy = True
                continue
            # Now copy the data to the BBP acc file
            pieces = line.split()
            if len(pieces) != 2:
                # We expect 2 values per line (time and acc)
                continue
            pieces = [float(piece) for piece in pieces]
            out_acc.write("%15.6e%15.6e%15.6e%15.6e\n" %
                          (pieces[0], pieces[1], pieces[1], pieces[1]))
        in_seis.close()
        out_acc.close()

        # Now create the vel bbp file!
        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(tmpdir,
                              "%d.%s.acc.000" % (self.sim_id, stat))
        ewfile = os.path.join(tmpdir,
                              "%d.%s.acc.090" % (self.sim_id, stat))
        udfile = os.path.join(tmpdir,
                              "%d.%s.acc.ver" % (self.sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s/wcc2bbp " % (self.install.A_GP_BIN_DIR) +
               "wcc2bbp=0 nsfile=%s ewfile=%s udfile=%s < %s >> %s 2>&1" %
               (nsfile, ewfile, udfile, out_bbp_acc_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
        # Now we need to integrate to get to velocity
        for comp in ['000', '090', 'ver']:
            file_in = os.path.join(tmpdir,
                                   "%d.%s.acc.%s" % (self.sim_id, stat, comp))
            file_out = os.path.join(tmpdir,
                                    "%s.%s.vel.%s" % (self.sim_id, stat, comp))
            bband_utils.check_path_lengths([file_in, file_out],
                                           bband_utils.GP_MAX_FILENAME)
            cmd = ("%s/integ_diff integ=1 filein=%s fileout=%s" %
                   (self.install.A_GP_BIN_DIR, file_in, file_out))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
        # Now we put together the file again as a velocity bbp file
        nsfile = os.path.join(tmpdir,
                              "%d.%s.vel.000" % (self.sim_id, stat))
        ewfile = os.path.join(tmpdir,
                              "%d.%s.vel.090" % (self.sim_id, stat))
        udfile = os.path.join(tmpdir,
                              "%d.%s.vel.ver" % (self.sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s/wcc2bbp " % (self.install.A_GP_BIN_DIR) +
               "wcc2bbp=1 nsfile=%s ewfile=%s udfile=%s > %s 2>> %s" %
               (nsfile, ewfile, udfile, out_bbp_vel_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

    def find_empirical_file(self):
        """
        This function finds the correct empirical file to use based on
        a file with magnitude ranges
        """
        mag = round(self.config.CFGDICT['magnitude'], 2)
        emp_amp_file = None

        infile = open(self.empirical_ranges, 'r')
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                # Skip comments
                continue
            pieces = line.split()
            m_min = float(pieces[0])
            m_max = float(pieces[1])
            # Check if event magnitude is within range
            if mag >= m_min and mag <= m_max:
                # Found the file we need
                emp_amp_file = pieces[2]
                break
        infile.close()

        if emp_amp_file is None:
            raise bband_utils.ParameterError("Cannot find empirical_amp file "
                                             "for event of magnitude %f" %
                                             (mag))
        return emp_amp_file

    def create_exsim_param_file(self):
        """
        This function creates the parameter file needed by ExSim
        """
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir_mod = os.path.join(self.install.A_TMP_DATA_DIR,
                                    str(self.sim_id),
                                    "exsim_%s" % (sta_base))
        a_param_template = os.path.join(a_indir, self.r_param_template)
        a_param_out = os.path.join(a_tmpdir_mod, self.config.PARAM_FILE)
        output_stem = "exsim-output-%s" % (sta_base)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        vmodel_params = vel_obj.get_codebase_params('exsim')

        # Look for additional files needed by ExSim
        if 'CRUSTAL_AMP' in vmodel_params:
            self.crustal_amp = os.path.join(vel_obj.base_dir,
                                            vmodel_params['CRUSTAL_AMP'])
        else:
            raise bband_utils.ParameterError("Cannot find crustal_amp "
                                             "parameter in velocity "
                                             "model %s" %
                                             (self.vmodel_name))
        if 'SITE_AMP' in vmodel_params:
            self.site_amp = os.path.join(vel_obj.base_dir,
                                         vmodel_params['SITE_AMP'])
        else:
            raise bband_utils.ParameterError("Cannot find site_amp "
                                             "parameter in velocity "
                                             "model %s" %
                                             (self.vmodel_name))
        if 'EMPIRICAL_DIR' in vmodel_params:
            self.empirical_dir = os.path.join(vel_obj.base_dir,
                                              vmodel_params['EMPIRICAL_DIR'])
        else:
            raise bband_utils.ParameterError("Cannot find empirical_dir "
                                             "parameter in velocity "
                                             "model %s" %
                                             (self.vmodel_name))
        if 'EMPIRICAL_RANGES' in vmodel_params:
            self.empirical_ranges = os.path.join(vel_obj.base_dir,
                                                 vmodel_params['EMPIRICAL_RANGES'])
        else:
            raise bband_utils.ParameterError("Cannot find empirical_ranges "
                                             "parameter in velocity "
                                             "model %s" %
                                             (self.vmodel_name))
        if 'SLIP_WEIGHTS' in vmodel_params:
            self.slip_weights = os.path.join(vel_obj.base_dir,
                                             vmodel_params['SLIP_WEIGHTS'])
        else:
            raise bband_utils.ParameterError("Cannot find slip_weights "
                                             "parameter in velocity "
                                             "model %s" %
                                             (self.vmodel_name))
        # Look for KAPPA and STRESS
        if 'KAPPA' in vmodel_params:
            self.kappa = float(vmodel_params['KAPPA'])
        if 'STRESS' in vmodel_params:
            self.stress = float(vmodel_params['STRESS'])

        # Check if we need to calculate stress
        if 'CALCULATE_STRESS' in vmodel_params:
            if float(vmodel_params['CALCULATE_STRESS']) == True:
                # Calculate stress based on depth of hypocenter
                self.stress = self.config.calculate_stress()

        # Stage these files into tmpdir_mod directory
        shutil.copy2(self.crustal_amp, a_tmpdir_mod)
        shutil.copy2(self.site_amp, a_tmpdir_mod)
        shutil.copy2(self.slip_weights, a_tmpdir_mod)

        # Now we need to figure out which empirical_amps file to use
        empirical_file = self.find_empirical_file()
        empirical_in = os.path.join(self.empirical_dir, empirical_file)
        empirical_out = os.path.join(a_tmpdir_mod, self.config.EMPIRICAL_AMPS)
        shutil.copy2(empirical_in, empirical_out)

        # Ok, need to write ExSim's param file now!
        self.template_in = open(a_param_template, 'r')
        self.param_out = open(a_param_out, 'w')

        # Replace parameters in the order they appear in ExSim's template
        self.template_replace(["<MAG>", "<STRESS>"],
                              [self.config.CFGDICT['magnitude'], self.stress])
        self.template_replace(["<MAG>", "<STRESS>", "<KAPPA>"],
                              [self.config.CFGDICT['magnitude'],
                               self.stress, self.kappa])
        self.template_replace(["<LAT>", "<LON>"],
                              [self.config.CFGDICT['lat_top_edge'],
                               self.config.CFGDICT['lon_top_edge']])
        self.template_replace(["<STRIKE>", "<DIP>", "<DEPTH>"],
                              [self.config.CFGDICT['strike'],
                               self.config.CFGDICT['dip'],
                               self.config.CFGDICT['depth_to_top']])
        self.template_replace(["<F_LEN>", "<F_WID>"],
                              [self.config.CFGDICT["fault_length"],
                               self.config.CFGDICT["fault_width"]])
        self.template_replace(["<HYPO_ALONG_STK>", "<HYPO_DOWN_DIP>"],
                              [self.config.CFGDICT["hypo_along_stk"],
                               self.config.CFGDICT["hypo_down_dip"]])
        self.template_replace(["<OUTPUT_STEM>"], [output_stem])
        self.template_replace(["<CRUSTAL_AMP_FILE>"],
                              [os.path.basename(self.crustal_amp)])
        self.template_replace(["<SITE_AMP_FILE>"],
                              [os.path.basename(self.site_amp)])
        self.template_replace(["<EMPIRICAL_AMP_FILE>"],
                              [os.path.basename(self.config.EMPIRICAL_AMPS)])
        self.template_replace(["<SEED>"],
                              [int(self.config.CFGDICT['seed'])])
        self.template_replace(["<SLIP_WEIGHTS_FILE>"],
                              [os.path.basename(self.slip_weights)])
        self.template_replace(["<NUM_STATIONS>"], [self.num_stations])
        # Now, copy everything else!
        self.template_replace(["<END_OF_TEMPLATE_FILE>"], [""])

        # Almost done, now we need to add stations
        site_list = self.stat_list.get_station_list()

        # Check for maximum number of stations
        if len(site_list) > self.config.MAX_STATIONS:
            raise bband_utils.ParameterError("Too many stations in "
                                             "the station list: %d. " %
                                             (len(site_list)) +
                                             "Maximum limit is %d." %
                                             (self.config.MAX_STATIONS))

        for site in site_list:
            self.param_out.write("%f  %f\n" % (site.lat, site.lon))

        # Done! Close files.
        self.template_in.close()
        self.param_out.close()

        # Make copy of ExSIM param file in input directory
        shutil.copy2(a_param_out, a_indir)

    def run(self):
        """
        This function prepares the parameter file for ExSim, invokes
        it, and formats its output to be compatible with the Broadband
        Platform
        """
        print("ExSIM".center(80, '-'))

        self.install = InstallCfg.getInstance()
        self.config = ExSimCfg(os.path.join(self.install.A_IN_DATA_DIR,
                                            str(self.sim_id),
                                            self.r_srcfile))
        install = self.install
        config = self.config
        sim_id = self.sim_id

        self.kappa = config.KAPPA
        self.stress = config.STRESS
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.exsim_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR,
                                    str(sim_id),
                                    "exsim_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")
        #
        # Make sure the output and two tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir, a_param_outdir],
                           print_cmd=False)

        a_stations = os.path.join(a_indir, self.r_stations)
        self.stat_list = StationList(a_stations)
        self.num_stations = len(self.stat_list.site_list)

        # Need to create ExSim's parameter file
        self.create_exsim_param_file()

        # Keep a copy of ExSim's parameter file in outdata
        a_param_file = os.path.join(a_tmpdir_mod, self.config.PARAM_FILE)
        shutil.copy2(a_param_file, os.path.join(a_param_outdir,
                                                self.config.PARAM_FILE))

        # Run in tmpdir subdir to isolate temp fortran files
        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        exsim_bin = os.path.join(install.A_UWO_BIN_DIR,
                                 "exsim14")
        cmd = ("%s >> %s 2>&1 << END\n" % (exsim_bin, self.log) +
               "%s\n" % (self.config.PARAM_FILE) +
               "END")
        print("==> Running ExSIM...")
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Restore working directory
        os.chdir(old_cwd)

        # Need to copy and re-format output files
        output_stem = "exsim-output-%s" % (sta_base)
        site_list = self.stat_list.get_station_list()

        for site, idx in zip(site_list, range(1, self.num_stations + 1)):
            stat = site.scode
            #out_rotd50_base = "%d.%s.rd50" % (sim_id, stat)
            #out_rotd50 = os.path.join(a_outdir, out_rotd50_base)
            #in_exsim_base = "%s_psa_fs_s%03d.out" % (output_stem, idx)
            #in_exsim_file = os.path.join(a_tmpdir_mod, "PSA", in_exsim_base)
            #print("==> Writing RotD50 file for station: %s" % (stat))
            #self.convert_exsim_to_rd50(in_exsim_file, out_rotd50)
            in_exsim_acc_base = "%sS%03diter001.acc" % (output_stem, idx)
            in_exsim_acc_file = os.path.join(a_tmpdir_mod,
                                             "ACC",
                                             in_exsim_acc_base)
            print("==> Writing BBP file for station: %s" % (stat))
            self.convert_exsim_to_bbp(in_exsim_acc_file, stat,
                                      a_tmpdir, a_outdir)

        print("ExSIM Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = ExSim(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
               sim_id=int(sys.argv[5]))
    ME.run()
