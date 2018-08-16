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

Broadband Platform Version of Rob Graves jbsim script
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import shutil

# Import Broadband modules
import bband_utils
import velocity_models
from match_cfg import MatchCfg
from install_cfg import InstallCfg
from station_list import StationList

class Match(object):
    """
    Implement Rob Graves match.csh as a python component
    """

    def __init__(self, i_r_stations, i_vmodel_name, sim_id=0,
                 acc=True, pow2=False):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.vmodel_name = i_vmodel_name
        self.acc = acc
        self.pow2 = pow2
        self.phase = None
        self.hf_fhi = None
        self.lf_flo = None

    def run(self):
        """
        Runs the match module to merge low and high frequency seismograms
        """
        print("Match".center(80, '-'))

        install = InstallCfg.getInstance()
        config = MatchCfg()

        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.match_%s.log" % (sim_id, sta_base))

        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        # Make sure tmpdir exists
        dirs = [a_tmpdir]
        bband_utils.mkdirs(dirs, print_cmd=False)

        pow2_param = 0
        if self.pow2:
            pow2_param = 1

        # Start with defaults
        self.phase = config.PHASE
        self.hf_fhi = config.HF_FHI
        self.lf_flo = config.LF_FLO

        # Set match method
        if config.MATCH_METHOD == 1:
            self.phase = 1
        elif config.MATCH_METHOD == 2:
            val = 1.0 / (2.0 * config.HF_ORD)
            self.hf_fhi = (self.hf_fhi *
                           math.exp(val * math.log(math.sqrt(2.0) - 1.0)))
            val = -1.0 / (2.0 * config.LF_ORD)
            self.lf_flo = (self.lf_flo *
                           math.exp(val * math.log(math.sqrt(2.0) - 1.0)))

        #
        # Read and parse the station list with this call
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))

        # Check for velocity model-specific parameters
        vmodel_params = vel_obj.get_codebase_params('gp')

        # Figure out what DT we should use when resampling

        # Figure out the LF DT value
        if self.acc:
            seis_ext = '.acc.bbp'
        else:
            seis_ext = '.bbp'

        lf_seis = None
        hf_seis = None

        # Find one LF seismogram and one HF seismogram
        for sites in site_list:
            site = sites.scode
            if os.path.exists(os.path.join(a_tmpdir,
                                           "%d.%s-lf%s" %
                                           (sim_id, site,
                                            seis_ext))):
                lf_seis = os.path.join(a_tmpdir,
                                       "%d.%s-lf%s" %
                                       (sim_id, site,
                                        seis_ext))
                if os.path.exists(os.path.join(a_tmpdir,
                                               "%d.%s-hf%s" %
                                               (sim_id, site,
                                                seis_ext))):
                    hf_seis = os.path.join(a_tmpdir,
                                           "%d.%s-hf%s" %
                                           (sim_id, site,
                                            seis_ext))
                    break

        # Need one of each
        if lf_seis is None:
            raise bband_utils.ParameterError("Cannot find a LF seismogram")
        if hf_seis is None:
            raise bband_utils.ParameterError("Cannot find a HF seismogram")

        # Pick DT from these files
        lf_dt = None
        lf_file = open(lf_seis)
        for line in lf_file:
            line = line.strip()
            if line.startswith("#") or line.startswith("%"):
                continue
            # Got to first timestamp. Now, pick two consecutive
            # timestamps values
            lf_t1 = float(line.strip().split()[0])
            lf_t2 = float(lf_file.next().strip().split()[0])
            # Subtract the two times
            lf_dt = lf_t2 - lf_t1
            # All done!
            break
        lf_file.close()

        if lf_dt is None:
            raise bband_utils.ParameterError("Cannot find LF_DT!")

        hf_dt = None
        hf_file = open(hf_seis)
        for line in hf_file:
            line = line.strip()
            if line.startswith("#") or line.startswith("%"):
                continue
            # Got to first timestamp. Now, pick two consecutive
            # timestamps values
            hf_t1 = float(line.strip().split()[0])
            hf_t2 = float(hf_file.next().strip().split()[0])
            # Subtract the two times
            hf_dt = hf_t2 - hf_t1
            # All done!
            break
        hf_file.close()

        if hf_dt is None:
            raise bband_utils.ParameterError("Cannot find HF_DT!")

        # In the GP method, we can potentially have two independent DT
        # values, one used by the rupture generator and the
        # low-frequency jbsim seismogram simulator, and another value
        # used by the high-frequency hfsims program. We have to use
        # the smaller of these two values in order to properly combine
        # the low-, and high-frequency seismograms.

        new_dt = min(lf_dt, hf_dt)

        # Go through the stations
        for sites in site_list:
            # Pick station name
            site = sites.scode
            #
            # We have a verbose of silent invocation. This is a very
            # verbose program so our default writes to dev/null
            #

            #
            # There are multiple possibilities; either we have
            # separate HF and LF files, we have HF and .bbp, LF and
            # .bbp, or just .bbp.  In all cases, we need to separate
            # them to get components.
            #
            hf_exists = False
            lf_exists = False

            if not self.acc:
                print("==> Processing velocity seismograms for station: %s" %
                      (site))
                # Need to convert to acc first
                if os.path.exists(os.path.join(a_tmpdir,
                                               "%d.%s-hf.bbp" %
                                               (sim_id, site))):
                    hf_exists = True
                if os.path.exists(os.path.join(a_tmpdir,
                                               "%d.%s-lf.bbp" %
                                               (sim_id, site))):
                    lf_exists = True

                # If no files exist for this station, make a note and continue
                if not hf_exists and not lf_exists:
                    print("===> No velocity seismograms found!")
                    print("===> Skipping station...")
                    continue

                # First process HF files to convert velocity to acceleration

                # Create path names and check if their sizes are
                # within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.000" % (sim_id, site))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.090" % (sim_id, site))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.ver" % (sim_id, site))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s-hf.bbp" % (sim_id, site))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                # Run wcc2bbp
                cmd = ("%s " %
                       (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "wcc2bbp=0 < %s >> %s 2>&1" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                for comp in config.COMPS:
                    # Create path names and check if their sizes
                    # are within bounds
                    filein = os.path.join(a_tmpdir,
                                          "%d.%s-hf.%s" %
                                          (sim_id, site, comp))
                    fileout = os.path.join(a_tmpdir,
                                           "%d.%s-hf.acc.%s" %
                                           (sim_id, site, comp))

                    bband_utils.check_path_lengths([filein, fileout],
                                                   bband_utils.GP_MAX_FILENAME)

                    cmd = ("%s diff=1 " %
                           (os.path.join(install.A_GP_BIN_DIR,
                                         "integ_diff")) +
                           "filein=%s fileout=%s" % (filein, fileout))
                    bband_utils.runprog(cmd, abort_on_error=True,
                                        print_cmd=False)

                # Create path names and check if their sizes are within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.acc.000" % (sim_id, site))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.acc.090" % (sim_id, site))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.acc.ver" % (sim_id, site))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s-hf.acc.bbp" % (sim_id, site))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s " %
                       (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                # Then process LF files to convert velocity to acceleration

                # Create path names and check if their sizes are within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.000" % (sim_id, site))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.090" % (sim_id, site))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.ver" % (sim_id, site))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s-lf.bbp" % (sim_id, site))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s " %
                       (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "wcc2bbp=0 < %s >> %s 2>&1" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                for comp in config.COMPS:
                    # Create path names and check if their sizes
                    # are within bounds
                    filein = os.path.join(a_tmpdir,
                                          "%d.%s-lf.%s" %
                                          (sim_id, site, comp))
                    fileout = os.path.join(a_tmpdir,
                                           "%d.%s-lf.acc.%s" %
                                           (sim_id, site, comp))

                    bband_utils.check_path_lengths([filein, fileout],
                                                   bband_utils.GP_MAX_FILENAME)

                    cmd = ("%s " %
                           (os.path.join(install.A_GP_BIN_DIR,
                                         "integ_diff")) +
                           "diff=1 filein=%s fileout=%s" %
                           (filein, fileout))
                    bband_utils.runprog(cmd, abort_on_error=True,
                                        print_cmd=False)

                # Create path names and check if their sizes are within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.acc.000" % (sim_id, site))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.acc.090" % (sim_id, site))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.acc.ver" % (sim_id, site))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s-lf.acc.bbp" % (sim_id, site))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s " %
                       (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # We should have acceleration files at this point
            hf_exists = False
            lf_exists = False

            if os.path.exists(os.path.join(a_tmpdir,
                                           "%d.%s-hf.acc.bbp" %
                                           (sim_id, site))):
                hf_exists = True
            if os.path.exists(os.path.join(a_tmpdir,
                                           "%d.%s-lf.acc.bbp" %
                                           (sim_id, site))):
                lf_exists = True

            print("==> Processing acceleration seismograms for station: %s" %
                  (site))

            # If no files exist for this station, make a note and continue
            if not hf_exists and not lf_exists:
                print("===> No acceleration seismograms found!")
                print("===> Skipping station...")
                continue

            #
            # Convert HF file to wcc components
            #

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s-hf.acc.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s-hf.acc.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s-hf.acc.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s-hf.acc.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          "nsfile=%s ewfile=%s udfile=%s " %
                          (nsfile, ewfile, udfile) +
                          "wcc2bbp=0 < %s >> %s 2>&1" %
                          (bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

            #
            # Convert LF file to wcc components
            #

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s-lf.acc.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s-lf.acc.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s-lf.acc.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s-lf.acc.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          "nsfile=%s ewfile=%s udfile=%s " %
                          (nsfile, ewfile, udfile) +
                          "wcc2bbp=0 < %s >> %s 2>&1" %
                          (bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

            #
            # Process each component
            #
            for entries in config.COMPS:
                compo = entries

                #
                # HF First
                #
                listfile = os.path.join(a_tmpdir, "%s.%s.hf.%s" %
                                        (config.FILTLIST, sta_base, compo))
                bband_utils.check_path_lengths([listfile],
                                               bband_utils.GP_MAX_FILENAME)

                # Create wcc_tfilter input file
                out = open(listfile, 'w')
                # Contains HF input file
                infile = os.path.join(a_tmpdir,
                                      "%d.%s-hf.acc.%s" %
                                      (sim_id, site, compo))
                out.write("%s\n" % infile)
                out.flush()
                out.close()

                # Also check infile
                bband_utils.check_path_lengths([infile],
                                               bband_utils.GP_MAX_FILENAME)

                #
                # Pre-filter and resample HF file
                #
                shutil.copy2(infile, "%s.prefilter" % infile)
                progstring = ("%s " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_tfilter")) +
                              "filelist=%s order=%d fhi=%f flo=%s " %
                              (listfile, config.HF_ORD, self.hf_fhi,
                               config.HF_FLO) +
                              "inbin=0 outbin=0 phase=%d " %
                              (self.phase) +
                              "outpath=%s >> %s 2>&1" %
                              (a_tmpdir, self.log))
                bband_utils.runprog(progstring, abort_on_error=True,
                                    print_cmd=False)

                outfile = os.path.join(a_tmpdir, "%d.%s-hf-resamp.%s" %
                                       (sim_id, site, compo))
                bband_utils.check_path_lengths([outfile],
                                               bband_utils.GP_MAX_FILENAME)

                progstring = ("%s newdt=%f " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_resamp_arbdt"), new_dt) +
                              "pow2=%d infile=%s outfile=%s >> %s 2>&1" %
                              (pow2_param, infile, outfile, self.log))
                bband_utils.runprog(progstring, abort_on_error=True,
                                    print_cmd=False)

                #
                # LF Next
                #
                listfile = os.path.join(a_tmpdir, "%s.%s.lf.%s" %
                                        (config.FILTLIST, sta_base, compo))
                bband_utils.check_path_lengths([listfile],
                                               bband_utils.GP_MAX_FILENAME)

                # Create wcc_tfilter input file
                out = open(listfile, 'w')
                # Contains LF input file
                infile = os.path.join(a_tmpdir,
                                      "%d.%s-lf.acc.%s" %
                                      (sim_id, site, compo))
                out.write("%s\n" % infile)
                out.flush()
                out.close()

                # Also check infile
                bband_utils.check_path_lengths([infile],
                                               bband_utils.GP_MAX_FILENAME)

                #
                # Pre-filter and resample LF file
                #
                shutil.copy2(infile, "%s.prefilter" % infile)
                progstring = ("%s " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_tfilter")) +
                              "filelist=%s order=%d fhi=%f flo=%s " %
                              (listfile, config.LF_ORD, config.LF_FHI,
                               self.lf_flo) +
                              "inbin=0 outbin=0 phase=%d " %
                              (self.phase) +
                              "outpath=%s >> %s 2>&1 " %
                              (a_tmpdir, self.log))
                bband_utils.runprog(progstring, print_cmd=False)

                outfile = os.path.join(a_tmpdir, "%d.%s-lf-resamp.%s" %
                                       (sim_id, site, compo))
                bband_utils.check_path_lengths([outfile],
                                               bband_utils.GP_MAX_FILENAME)

                progstring = ("%s " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_resamp_arbdt")) +
                              "newdt=%f pow2=%d " %
                              (new_dt, pow2_param) +
                              "infile=%s outfile=%s >> %s 2>&1" %
                              (infile, outfile, self.log))
                bband_utils.runprog(progstring, abort_on_error=True,
                                    print_cmd=False)

                #
                # Add LF and HF resampled acc seismograms
                #

                # Check all path lengths
                infile1 = os.path.join(a_tmpdir, "%d.%s-lf-resamp.%s" %
                                       (sim_id, site, compo))
                infile2 = os.path.join(a_tmpdir, "%d.%s-hf-resamp.%s" %
                                       (sim_id, site, compo))
                outfile = os.path.join(a_tmpdir, "%d.%s.acc.add.%s" %
                                       (sim_id, site, compo))
                bband_utils.check_path_lengths([infile1, infile2, outfile],
                                               bband_utils.GP_MAX_FILENAME)

                progstring = ("%s " %
                              (os.path.join(install.A_GP_BIN_DIR, "wcc_add")) +
                              "f1=1.00 t1=%f inbin1=0 infile1=%s " %
                              (config.LF_TSTART, infile1) +
                              "f2=1.00 t2=%f inbin2=0 infile2=%s " %
                              (config.HF_TSTART, infile2) +
                              "outbin=0 outfile=%s >> %s 2>&1" %
                              (outfile, self.log))
                bband_utils.runprog(progstring, abort_on_error=True,
                                    print_cmd=False)

                #
                # Create combined velocity files
                #

                # Check path lengths
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.acc.add.%s" %
                                      (sim_id, site, compo))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.%s" %
                                       (sim_id, site, compo))
                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s integ=1 filein=%s fileout=%s" %
                       (os.path.join(install.A_GP_BIN_DIR,
                                     "integ_diff"), filein, fileout))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # We have all the component files, create velocity seismogram

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s wcc2bbp=1 " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          'title="Sim NGAH, stat=%s" ' % site +
                          'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                          (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

            # Copy velocity bbp file to outdir
            shutil.copy2(os.path.join(a_tmpdir, "%d.%s.bbp" %
                                      (sim_id, site)),
                         os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                      (sim_id, site)))

            # Also create acceleration bbp file in outdir

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "wcc2bbp=0 < %s >> %s 2>&1" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            for comp in config.COMPS:
                # Create path names and check if their sizes are within bounds
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.%s" %
                                      (sim_id, site, comp))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.acc.%s" %
                                       (sim_id, site, comp))

                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s diff=1 filein=%s fileout=%s" %
                       (os.path.join(install.A_GP_BIN_DIR,
                                     "integ_diff"), filein, fileout))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.acc.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s " % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Copy acceleration bbp file to outdir
            shutil.copy2(os.path.join(a_tmpdir, "%d.%s.acc.bbp" %
                                      (sim_id, site)),
                         os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                      (sim_id, site)))

        print("Match Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = Match(sys.argv[1], sys.argv[2], int(sys.argv[3]), acc=False)
    ME.run()
