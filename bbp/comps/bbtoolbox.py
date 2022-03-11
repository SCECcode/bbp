#!/bin/env python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Broadband Platform Version of Martin Mai BBcoda2.csh
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import random
import shutil

# Import Broadband modules
import bband_utils
import stas2files
import velocity_models
from station_list import StationList
from install_cfg import InstallCfg
from bbtoolbox_cfg import BBToolboxCfg
from bbtoolbox_correlation import generate_matrices

class BBToolbox(object):

    def __init__(self, i_r_scattering, i_r_velmodel, i_r_srcfile,
                 i_r_srffile, i_r_stations, vmodel_name, sim_id=0,
                 **kwargs):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srcfiles = []
        self.r_scattering = i_r_scattering
        self.r_srffile = i_r_srffile
        self.r_xyz_srffile = 'xyz_' + i_r_srffile
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.config = None
        self.iseed = None
        self.fmax = None
        self.kappa = None
        self.q_coda = None
        self.fdec = None
        self.source_func = None
        self.gs_flag = None
        self.ngaw_flag = None
        self.tr_sca = None
        self.afac = None
        self.bfac = None
        self.str_fac = None
        self.correlation_file = None
        self.corr_flag = None

        # Get all src files that were passed to us
        if kwargs is not None and len(kwargs) > 0:
            for idx in range(len(kwargs)):
                self.r_srcfiles.append(kwargs['src%d' % (idx)])
        else:
            # Not a multisegment run, just use the single src file
            self.r_srcfiles.append(i_r_srcfile)

    def create_bbtoolbox_files(self, stat_file):
        """
        This function creates the files needed by bbtoolbox, including
        the scattering file (if not provided), the station file, and
        the parameter file
        """
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(self.sim_id))
        a_tmpdir_mod = os.path.join(self.install.A_TMP_DATA_DIR,
                                    str(self.sim_id),
                                    "bbtoolbox_%s" % (sta_base))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")

        stat_list = StationList(stat_file)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        vmodel_params = vel_obj.get_codebase_params('sdsu')

        # Look for the source function parameter
        if 'SOURCE_FUNC' in vmodel_params:
            self.source_func = vmodel_params['SOURCE_FUNC']

        # Set up correlation parameter
        self.corr_flag = self.config.corr_flag
        # Create correlation matrices
        if self.corr_flag > 0:
            generate_matrices(self.install.A_SDSU_DATA_DIR, a_tmpdir,
                              stat_file, a_tmpdir_mod, self.sim_id)

        # Take care of scattering file
        if not self.r_scattering:
            # Need to create our file
            scattering_template = os.path.join(self.install.A_SDSU_DATA_DIR,
                                               "scattering_generic.dat")
            self.r_scattering = "scattering.dat"
            a_scattering = os.path.join(a_indir, self.r_scattering)

            # Look for KAPPA
            if 'KAPPA' in vmodel_params:
                self.kappa = float(vmodel_params['KAPPA'])
            # Look for FMAX
            if 'FMAX' in vmodel_params:
                self.fmax = float(vmodel_params['FMAX'])
            if 'Q' in vmodel_params:
                self.q_coda = float(vmodel_params['Q'])
            if 'FDEC' in vmodel_params:
                self.fdec = float(vmodel_params['FDEC'])
            if 'GS_FLAG' in vmodel_params:
                self.gs_flag = float(vmodel_params['GS_FLAG'])
            if 'NGAW_FLAG' in vmodel_params:
                self.ngaw_flag = float(vmodel_params['NGAW_FLAG'])
            if 'TR_SCA' in vmodel_params:
                self.tr_sca = float(vmodel_params['TR_SCA'])
            if 'AFAC' in vmodel_params:
                self.afac = float(vmodel_params['AFAC'])
            if 'BFAC' in vmodel_params:
                self.bfac = float(vmodel_params['BFAC'])
            if 'STR_FAC' in vmodel_params:
                self.str_fac = float(vmodel_params['STR_FAC'])

            # Check if we need to calculate stress
            if 'CALCULATE_STRESS' in vmodel_params:
                if float(vmodel_params['CALCULATE_STRESS']) == True:
                    # Calculate stress based on depth of hypocenter
                    self.str_fac = self.config.calculate_stress()

            # Open template and output files
            scat_in = open(scattering_template, 'r')
            scat_out = open(a_scattering, 'w')

            for line in scat_in:
                if line.find(r"\* iseed - seed number for scattering") >= 0:
                    # This is the iseed line, insert the random iseed here
                    pos = line.find(r"\* iseed - seed number for scattering")
                    scat_out.write("%d   %s" %
                                   (self.iseed,
                                    line[pos:]))
                elif line.find(r"\* kappa - kappa at the site") >= 0:
                    # This is the kappa line, insert self.kappa here
                    pos = line.find(r"\* kappa - kappa at the site")
                    scat_out.write("%.3f   %s" %
                                   (self.kappa,
                                    line[pos:]))
                elif line.find(r"\* fmax - ") >= 0:
                    # This is the fmax line, insert self.fmax here
                    pos = line.find(r"\* fmax - ")
                    scat_out.write("%.2f   %s" %
                                   (self.fmax,
                                    line[pos:]))
                elif line.find(r"\* Q - Q for the coda") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* Q - Q for the coda")
                    scat_out.write("%.1f   %s" %
                                   (self.q_coda,
                                    line[pos:]))
                elif line.find(r"\* fdec - see equation") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* fdec - see equation")
                    scat_out.write("%.2f    %s" %
                                   (self.fdec,
                                    line[pos:]))
                elif line.find(r"\* gs_flag - determine type") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* gs_flag - determine type")
                    scat_out.write("%d    %s" %
                                   (int(self.gs_flag),
                                    line[pos:]))
                elif line.find(r"\* ngaw_flag - GMPEs") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* ngaw_flag - GMPEs")
                    scat_out.write("%d    %s" %
                                   (int(self.ngaw_flag),
                                    line[pos:]))
                elif line.find(r"\* Tr_sca - scaling factor") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* Tr_sca - scaling factor")
                    scat_out.write("%.4f  %s" %
                                   (self.tr_sca,
                                    line[pos:]))
                elif line.find(r"\* afac - qk factor") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* afac - qk factor")
                    scat_out.write("%.1f    %s" %
                                   (self.afac,
                                    line[pos:]))
                elif line.find(r"\* bfac - qk factor") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* bfac - qk factor")
                    scat_out.write("%.1f    %s" %
                                   (self.bfac,
                                    line[pos:]))
                elif line.find(r"\* str_fac - Brune stress") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* str_fac - Brune stress")
                    scat_out.write("%.2e %s" %
                                   (self.str_fac,
                                    line[pos:]))
                elif line.find(r"\* cseed - seed number") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* cseed - seed number")
                    scat_out.write("%d   %s" %
                                   (self.config.SEED,
                                   line[pos:]))
                elif line.find(r"\* corr_flag") >= 0:
                    # This is the line, insert here
                    pos = line.find(r"\* corr_flag")
                    scat_out.write("%d    %s" %
                                   (int(self.corr_flag),
                                    line[pos:]))
                else:
                    scat_out.write(line)

            # Done
            scat_in.close()
            scat_out.flush()
            scat_out.close()

        # Keep copy of scattering file in outdata
        shutil.copy2(a_scattering, os.path.join(a_param_outdir,
                                                self.r_scattering))

        # Convert station file
        a_tmpfile = "station_%s.coords" % (sta_base)
        a_sdsu_stat_list = os.path.join(a_tmpdir_mod,
                                        "bbtstations_%s.tmp" % (sta_base))
        a_sdsu_extended_fault = os.path.join(a_indir, "extended_fault")
        param_filename = stas2files.bbp2sdsu_statlist(a_indir, stat_list,
                                                      a_sdsu_stat_list,
                                                      self.r_srffile,
                                                      self.r_xyz_srffile,
                                                      a_sdsu_extended_fault,
                                                      a_tmpfile)
        r_faultfile = os.path.basename(a_sdsu_extended_fault)
        # param_filename = stas2files.bbp2sdsu_statlist(a_indir, stat_list,
        #                                               a_sdsu_stat_list, hypo)
        # now a_sdsu_stat_list has X Y name vs rho kappa
        # a_sdsu_stat_list.par has bbextension, bbstat, bbhypo

        # Build real station list
        self.r_stations = "bbtstations_%s.dat" % (sta_base)
        stalist_fp = open(os.path.join(a_indir, self.r_stations), 'w')
        # write headers
        stalist_fp.write("/* STATIONS FILE FOR BROAD-BAND COMPUTATION CODE " +
                         "(P.M. MAI & K.B.OLSEN) */\n")
        stalist_fp.write("/* STATIONS COORDINATES ARE IN THE X-Y SYSTEM " +
                         "REPORTED IN FIG.1 OF APPENDIX A */\n\n")
        stalist_fp.write("/* INPUT DIRECTORY */\n")
        # Create input directory and file prefix for the stations files
        file_prefix = os.path.join(a_tmpdir_mod, "%d." % (self.sim_id))
        stalist_fp.write("%s\n\n" % (file_prefix))
        stalist_fp.write("/* FILES FORMAT [RGF BIN CMP 3SF] */\n")
        stalist_fp.write("\t3SF\n\n")
        stalist_fp.write("/* FILES EXTENSION OR BINARY FILE NAME */\n")
        glob_stat = "%s/*-lf.bbp" % (a_tmpdir)
        bbp_list = glob.glob(glob_stat)
        # Now, figure out the file suffix
        if len(bbp_list) > 0:
            file_suffix = "-lf.bbp"
        else:
            file_suffix = ".bbp"

        # Write suffix
        stalist_fp.write("%s\n\n" % (file_suffix))

        # Write header for station list
        stalist_fp.write("/*\tX\tY\tNAME\tVs\tRho\tKappa */\n")

        # Now, append the station list we have in a_sdsu_stat_list
        conv_list_fp = open(a_sdsu_stat_list, 'r')
        for line in conv_list_fp:
            stalist_fp.write(line)
            # Figure out if station file path is too long
            pieces = line.split()
            st_name = pieces[2]
            total_length = len(file_prefix) + len(st_name) + len(file_suffix)
            if total_length >= bband_utils.SDSU_MAX_FILENAME:
                # Close files
                stalist_fp.close()
                conv_list_fp.close()
                raise ValueError("station path for %s " % (st_name) +
                                 " is %d characters long, maximum is %d" %
                                 (total_length, bband_utils.SDSU_MAX_FILENAME))
        # Flush all data, and close this file
        stalist_fp.flush()
        stalist_fp.close()
        # Close station file
        conv_list_fp.close()

        # Keep copy of station file in outdata
        shutil.copy2(os.path.join(a_indir, self.r_stations),
                     os.path.join(a_param_outdir, self.r_stations))

        # Read param file
        conv_par_fp = open(param_filename, 'r')
        conv_par_data = conv_par_fp.readlines()
        conv_par_fp.close()

        # 2nd line is hypo coordinates
        hypo_line = conv_par_data[1].split(':')[1]
        hypo_coords = []
        for i in range(0, 3):
            hypo_coords.append(hypo_line.split()[i])
        min_box_dims = []
        min_box_line = conv_par_data[0].split(':')[1]
        for i in range(0, 2):
            min_box_dims.append(float(min_box_line.split()[i]))

        # FS: Feb-2013: Get magnitude directly from SRC file
        # FS: Mar-2013: We use this magnitude only when we don't have
        # a SRC file
        # get magnitude from 3rd line
        magnitude = float(conv_par_data[2].split(':')[1])

        self.r_bbparfile = "%d_%s.bbpar" % (self.sim_id, sta_base)
        parfile_name = os.path.join(a_indir, self.r_bbparfile)
        parfile_fp = open(parfile_name, 'w')
        parfile_fp.write("/* MODALITY FLAG: [0] LF-HF MERGING, " +
                         "[1] LF-SCATTERING, [2] LF-ISOCHRONE */\n")
        parfile_fp.write(" %d\n" % (self.config.MODALITY))
        parfile_fp.write("/* OUTPUT DIRECTORY */\n")
        parfile_fp.write('"%s"\n' % a_tmpdir_mod)
        parfile_fp.write('/* VELOCITY MODEL FILE (3D MODEL OR 1D MODEL) */\n')
        parfile_fp.write('"%s"\n' %
                         (os.path.join(a_indir, self.r_velmodel)))
        parfile_fp.write("/* STATIONS FILE REPORTING [X-Y] COORDINATES, " +
                         "FILENAMES AND PARAMETERS */\n")
        parfile_fp.write('"%s"\n' %
                         (os.path.join(a_indir, self.r_stations)))
        parfile_fp.write("/* OPTIONAL 2ND STATIONS FILE REPORTING ONLY " +
                         "FILENAMES - ONLY FOR MODALITY = 0  */\n")
        parfile_fp.write("2ndstations.dat\n")
        parfile_fp.write("/* FAULT MODEL TYPE: [POINT], " +
                         "[EXTENDED FAULT-MODEL FILE] */\n")
        parfile_fp.write(' extended "%s"\n' %
                         (os.path.join(a_indir, r_faultfile)))
        # parfile_fp.write(' point\n')
        parfile_fp.write("/* HYPOCENTER COORDINATES [X-Y-Z] IN KM */\n")
        parfile_fp.write("%.2f %.2f %.2f\n" % (float(hypo_coords[0]),
                                               float(hypo_coords[1]),
                                               float(hypo_coords[2])))
        parfile_fp.write('/* GRID DEFINITION [X-Y-Z] FOR RAYTRACING: ' +
                         '"NEAR-SIDE", GRID-SPACING (IN KM) */\n')
        parfile_fp.write("0.0 0.0 0.0 1.0\n")
        parfile_fp.write('/* GRID DEFINITION [X-Y-Z] FOR RAYTRACING: ' +
                         '"FAR-SIDE" (IN KM) */\n')
        if self.config.grid_x is not None and self.config.grid_y is not None:
            parfile_fp.write("%.1f %.1f %.1f\n" %
                             (self.config.grid_x,
                              self.config.grid_y,
                              self.config.grid_z))
        else:
            parfile_fp.write("%.1f %.1f %.1f\n" %
                             (round(min_box_dims[0] + 20.0, 0),
                              round(min_box_dims[1] + 20.0, 0),
                              self.config.grid_z))
        parfile_fp.write("/* SCATTERING PARAMETERS FILE */\n")
        parfile_fp.write('"%s"\n' %
                         (os.path.join(a_indir, self.r_scattering)))
        parfile_fp.write("/* EVENT MAGNITUDE */\n")
        if self.config.MAG is None:
            parfile_fp.write("%.2f\n" % (magnitude))
        else:
            parfile_fp.write("%.2f\n" % (self.config.MAG))
        parfile_fp.write("/* DOMINANT SOURCE MECHANISM [SS RS NS AL] */\n")
        parfile_fp.write("%s\n" % conv_par_data[3].split(":")[1].strip())
        parfile_fp.write("/* SOURCE TIME FUNCTION "
                         "[TRI BOX YOF DREG LIU NEW USER-DEF] */\n")
        parfile_fp.write("%s\n" % (self.source_func))
        parfile_fp.write("/* VERBOSE MODE [ON OFF] */\n")
        parfile_fp.write("off\n")
        parfile_fp.write("/* SRF FILE */\n")
        parfile_fp.write('"%s"\n' %
                         (os.path.join(a_indir, self.r_xyz_srffile)))
        parfile_fp.write("/* RAKE */\n")
        parfile_fp.write("%.2f\n" % (self.config.RAKE))
        parfile_fp.write("/* CORRELATION FILE FOR INTER-FREQUENCY CORRELATION*/\n")
        parfile_fp.write("Kinf.bin\n")
        parfile_fp.write("/* CORRELATION FILE FOR SPATIAL CORRELATION*/\n")
        parfile_fp.write("Ksp1.bin Ksp2.bin Ksp3.bin\n")
        parfile_fp.flush()
        parfile_fp.close()

        # Keep a copy in the outdata directory
        shutil.copy2(parfile_name, os.path.join(a_param_outdir,
                                                self.r_bbparfile))

    def run(self):
        """
        This function prepares the parameter file for BBToolbox, and
        then invokes it
        """
        print("SDSU BBToolBox".center(80, '-'))

        self.install = InstallCfg.getInstance()
        install = self.install
        sim_id = self.sim_id

        # Build path names
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.bbtoolbox_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR,
                                    str(sim_id),
                                    "bbtoolbox_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")
        #
        # Make sure the output and two tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir, a_param_outdir],
                           print_cmd=False)

        # Make sure BBToolbox works when starting from a srf file
        if self.r_srcfile:
            a_srcfile = os.path.join(a_indir, self.r_srcfile)
        else:
            a_srcfile = ""
        self.config = BBToolboxCfg(a_srcfile=a_srcfile)
        # Set default parameters
        config = self.config
        # Initialize random number with seed and calculate new iseed
        random.seed(config.SEED)
        self.iseed = int(random.random() * 10000)
        self.fmax = config.FMAX
        self.kappa = config.KAPPA
        self.q_coda = config.Q_CODA
        self.fdec = config.FDEC
        self.source_func = config.SOURCE_FUNC
        self.gs_flag = config.GS_FLAG
        self.ngaw_flag = config.NGAW_FLAG
        self.tr_sca = config.TR_SCA
        self.afac = config.AFAC
        self.bfac = config.BFAC
        self.str_fac = config.STR_FAC

        # Write valid par file, which includes correct references to
        # output dir, velocity model, stations list, fault
        # description, and scattering

        a_stations = os.path.join(a_indir, self.r_stations)

        # Need to create SDSU's BBToolbox files first
        self.create_bbtoolbox_files(a_stations)
        a_stations = os.path.join(a_indir,
                                  "bbtstations_%s.dat" %
                                  (sta_base))

        parfilename = os.path.join(a_tmpdir_mod, "parfilename_%s" % (sta_base))
        filename_fp = open(parfilename, "w")
        filename_fp.write('"%s"\n' % (os.path.join(a_indir, self.r_bbparfile)))
        filename_fp.flush()
        filename_fp.close()

        # Get list of stations
        stat_file_fp = open(a_stations, "r")
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

        # Check if we have stations in our list
        if len(stat_names) == 0:
            # No stations! We should output an error
            raise ValueError("No stations in the station list!")

        # Stagein seismograms to working dir
        for i in range(0, len(stat_names)):
            glob_list = glob.glob("%s/%d.%s*" % (a_tmpdir, sim_id,
                                                 stat_names[i]))
            for seismogram_file in glob_list:
                basename = os.path.basename(seismogram_file)
                shutil.copy2(seismogram_file, os.path.join(a_tmpdir_mod,
                                                           basename))

        # Run in tmpdir subdir to isolate temp fortran files
        os.chdir(a_tmpdir_mod)

        cmd = "%s/BBtoolbox.exe < %s >> %s 2>&1" % (install.A_SDSU_BIN_DIR,
                                                    parfilename, self.log)
        bband_utils.runprog(cmd, abort_on_error=True)

        for i in range(0, len(stat_names)):
            shutil.copy2("%s/BB.%s.hyb" % (a_tmpdir_mod, stat_names[i]),
                         "%s/%d.%s.bbp" % (a_tmpdir, sim_id, stat_names[i]))

        if config.copy_lf_seismograms:
            a_param_outdir = os.path.join(a_outdir, "param_files")
            for i in range(0, len(stat_names)):
                # Keep copy of lf seismogram files in outdata
                shutil.copy2("%s/%s.%s-lf.bbp" % (a_tmpdir, sim_id, stat_names[i]),
                             "%s/%s.%s-lf.bbp" % (a_param_outdir, sim_id, stat_names[i]))

        # Change to tmpdir so run.log file is put in tmpdir
        os.chdir(a_tmpdir)

        print("SDSU BBToolBox Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = BBToolbox(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                   sys.argv[5], None, sim_id=int(sys.argv[6]))
    ME.run()
