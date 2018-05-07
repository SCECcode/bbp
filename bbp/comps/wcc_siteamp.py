#!/usr/bin/env python
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

Broadband Platform Version of Rob Graves jbsim script
Outputs acceleration (cm/s/s)
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
from wcc_siteamp_cfg import WccSiteampCfg
from install_cfg import InstallCfg
from station_list import StationList

class WccSiteamp(object):
    """
    Implement Rob Graves siteamp.csh as a python component
    """

    def __init__(self, i_r_stations, method, i_vmodel_name, sim_id=0):
        """
        Initialize WccSiteamp parameters
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.method = method
        self.vmodel_name = i_vmodel_name
        self.config = None
        self.install = None
        self.log = None

    def process_separate_seismograms(self, site, sta_base, vs30, a_tmpdir):
        """
        Runs the site response module for separate low and high
        frequency seismograms
        """
        sim_id = self.sim_id
        config = self.config
        install = self.install

        flowcap = config.FLOWCAP

        # Run HF and LF components
        for freq in ['hf', 'lf']:
            print("**** Processing %s component..." % (freq))
            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "temp_%s.%s.000" % (sta_base, freq))
            ewfile = os.path.join(a_tmpdir,
                                  "temp_%s.%s.090" % (sta_base, freq))
            udfile = os.path.join(a_tmpdir,
                                  "temp_%s.%s.ver" % (sta_base, freq))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s-%s.bbp" % (sim_id, site,
                                                     freq))
            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            # Convert bbp files to 3 wcc files
            progstring = ("%s nsfile=%s ewfile=%s udfile=%s " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp"),
                           nsfile, ewfile, udfile) +
                          "wcc2bbp=0 < %s >> %s 2>&1" %
                          (bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

            for entries in config.COMPS:
                compo = entries

                # Create path names and check if their sizes
                # are within bounds
                filein = os.path.join(a_tmpdir,
                                      "temp_%s.%s.%s" %
                                      (sta_base, freq, compo))
                fileout = os.path.join(a_tmpdir,
                                       "temp_%s.acc.%s.%s" %
                                       (sta_base, freq, compo))
                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                # Convert to acceleration
                cmd = ("%s diff=1 " %
                       (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                       "filein=%s fileout=%s >> %s 2>&1" %
                       (filein, fileout, self.log))
                bband_utils.runprog(cmd, abort_on_error=True,
                                    print_cmd=False)

                # Get PGA for component from high frequency ONLY
                pga_filename = os.path.join(a_tmpdir,
                                            "%s.%s.pga.txt" % (site,
                                                               compo))

                # No need to check path lengths here
                progstring = ("%s < " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_getpeak")) +
                              "%s/temp_%s.acc.hf.%s > %s 2>> %s" %
                              (a_tmpdir, sta_base, compo,
                               pga_filename, self.log))
                # wcc_getpeak returns 33 even when it succeeds
                bband_utils.runprog(progstring, print_cmd=False)
                pga_file = open(pga_filename, "r")
                data = pga_file.readlines()
                pga_file.close()
                pga = float(data[0].split()[1]) / 981.0
                if freq == 'lf':
                    vref = config.LF_VREF
                else:
                    vref = config.HF_VREF

                # Create path names and check if their sizes
                # are within bounds
                filein = os.path.join(a_tmpdir,
                                      "temp_%s.acc.%s.%s" %
                                      (sta_base, freq, compo))
                fileout = os.path.join(a_tmpdir,
                                       "temp_%s.acc.amp.%s.%s" %
                                       (sta_base, freq, compo))
                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                # Pick the right model to use
                site_amp_model = config.SITEAMP_MODEL

                # Now, run site amplification
                progstring = ("%s pga=%f vref=%d " %
                              (os.path.join(install.A_GP_BIN_DIR,
                                            "wcc_siteamp14"), pga, vref) +
                              'vsite=%d model="%s" vpga=%d ' %
                              (vs30, site_amp_model,
                               config.HF_VREF) +
                              'flowcap=%f infile=%s outfile=%s ' %
                              (flowcap, filein, fileout) +
                              "fmidbot=%s fmin=%s >> %s 2>&1" %
                              (config.FMIDBOT, config.FMIN, self.log))
                bband_utils.runprog(progstring, abort_on_error=True)

                # Output becomes input
                filein = fileout
                fileout = os.path.join(a_tmpdir,
                                       "temp_%s.vel.amp.%s.%s" %
                                       (sta_base, freq, compo))
                bband_utils.check_path_lengths([fileout],
                                               bband_utils.GP_MAX_FILENAME)
                # Now integrate to get velocity
                cmd = ("%s integ=1 " %
                       (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                       "filein=%s fileout=%s >> %s 2>&1" %
                       (filein, fileout, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Do one time for acceleration

            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "temp_%s.acc.amp.%s.000" %
                                  (sta_base, freq))
            ewfile = os.path.join(a_tmpdir,
                                  "temp_%s.acc.amp.%s.090" %
                                  (sta_base, freq))
            udfile = os.path.join(a_tmpdir,
                                  "temp_%s.acc.amp.%s.ver" %
                                  (sta_base, freq))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s-%s.acc.bbp" %
                                   (sim_id, site, freq))
            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s wcc2bbp=1 units=cm/s/s " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          'title="%s Sim NGAH, stat=%s" ' %
                          (freq, site) +
                          'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                          (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)

            # Do it again for velocity

            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "temp_%s.vel.amp.%s.000" %
                                  (sta_base, freq))
            ewfile = os.path.join(a_tmpdir,
                                  "temp_%s.vel.amp.%s.090" %
                                  (sta_base, freq))
            udfile = os.path.join(a_tmpdir,
                                  "temp_%s.vel.amp.%s.ver" %
                                  (sta_base, freq))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s-%s.vel.bbp" %
                                   (sim_id, site, freq))
            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            progstring = ("%s wcc2bbp=1 units=cm/s " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          'title="%s Sim NGAH, stat=%s" ' %
                          (freq, site) +
                          'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                          (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)


    def process_hybrid_seismogram(self, site, sta_base, vs30,
                                  a_tmpdir, a_outdir):
        """
        Runs the site response module for a hybrid seismogram
        """
        sim_id = self.sim_id
        config = self.config
        install = self.install
        print("**** Reading broadband seismogram...")

        # Separate into HF and LF pieces:
        # 1) wcc2bbp
        # 2) convert to acceleration
        # 3) filter HF
        # 4) filter LF
        # 5) recombine HF and LF

        # Figure out vref to use in hybrid scenarios
        if self.method == "SDSU" or self.method == "UCSB":
            # For SDSU and UCSB we use whatever VREF was
            # used to calculate the LF seismograms
            vref = config.LF_VREF
        else:
            vref = config.GEN_ROCK_VS

        # Figure out where the input seismogram is located
        if self.method == "SDSU" or self.method == "UCSB":
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.bbp" % (sim_id, site))
        elif self.method == "EXSIM":
            bbpfile = os.path.join(a_outdir,
                                   "%d.%s.vel.bbp" %
                                   (sim_id, site))
        else:
            print("Unknown simulation method!")
            sys.exit(-1)

        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "temp_%s.000" % (sta_base))
        ewfile = os.path.join(a_tmpdir,
                              "temp_%s.090" % (sta_base))
        udfile = os.path.join(a_tmpdir,
                              "temp_%s.ver" % (sta_base))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        # Run wcc2bpp
        progstring = ("%s " %
                      (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                      "nsfile=%s " % (nsfile) +
                      "ewfile=%s " % (ewfile) +
                      "udfile=%s " % (udfile) +
                      "wcc2bbp=0 " +
                      "< %s >> %s 2>&1 " %
                      (bbpfile, self.log))

        bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

        for entries in config.COMPS:
            compo = entries
            # Create path names and check if their sizes are
            # within bounds
            filein = os.path.join(a_tmpdir,
                                  "temp_%s.%s" % (sta_base, compo))
            fileout = os.path.join(a_tmpdir,
                                   "temp_%s.acc.%s" % (sta_base, compo))
            bband_utils.check_path_lengths([filein, fileout],
                                           bband_utils.GP_MAX_FILENAME)

            # Convert from velocity to accl
            cmd = ("%s diff=1 " %
                   (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                   "filein=%s fileout=%s >> %s 2>&1" %
                   (filein, fileout, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            pga_filename = os.path.join(a_tmpdir,
                                        "%s.%s.pga.txt" % (site, compo))
            # No need to check path lengths here, it inputs
            # fileout from integ_diff
            filein = fileout
            progstring = ("%s < %s > %s 2>> %s" %
                          (os.path.join(install.A_GP_BIN_DIR,
                                        "wcc_getpeak"),
                           filein, pga_filename, self.log))
            # wcc_getpeak returns 33 even when it succeeds
            bband_utils.runprog(progstring, print_cmd=False)
            pga_file = open(pga_filename, "r")
            data = pga_file.readlines()
            pga_file.close()
            pga = float(data[0].split()[1]) / 981.0
            # Create path name for output file
            fileout = os.path.join(a_tmpdir,
                                   "temp_%s.amp.acc.%s" %
                                   (sta_base, compo))
            bband_utils.check_path_lengths([fileout],
                                           bband_utils.GP_MAX_FILENAME)

            # Pick the right model to use
            site_amp_model = config.SITEAMP_MODEL

            # Run site amp
            progstring = ("%s pga=%f vref=%d vsite=%d " %
                          (os.path.join(install.A_GP_BIN_DIR,
                                        "wcc_siteamp14"),
                           pga, vref, vs30) +
                          'model="%s" vpga=%d flowcap=%f ' %
                          (site_amp_model,
                           config.GEN_ROCK_VS, config.FLOWCAP) +
                          'infile=%s outfile=%s ' %
                          (filein, fileout) +
                          "fmidbot=%s fmin=%s >> %s 2>&1" %
                          (config.FMIDBOT, config.FMIN, self.log))
            bband_utils.runprog(progstring, abort_on_error=True)

            filein = fileout
            # Create path name for output file
            fileout = os.path.join(a_tmpdir,
                                   "temp_%s.vel.amp.%s" %
                                   (sta_base, compo))
            bband_utils.check_path_lengths([fileout],
                                           bband_utils.GP_MAX_FILENAME)
            # Convert to velocity
            cmd = ("%s integ=1 " %
                   (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                   "filein=%s fileout=%s >> %s 2>&1" %
                   (filein, fileout, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Create the velocity BBP file by combining the 3 vel components

        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "temp_%s.vel.amp.000" % (sta_base))
        ewfile = os.path.join(a_tmpdir,
                              "temp_%s.vel.amp.090" % (sta_base))
        udfile = os.path.join(a_tmpdir,
                              "temp_%s.vel.amp.ver" % (sta_base))
        bbpfile = os.path.join(a_tmpdir,
                               "%d.%s.vel.bbp" % (sim_id, site))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        progstring = ("%s wcc2bbp=1 units=cm/s " %
                      (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                      'title="Sim NGAH, stat=%s" ' % (site) +
                      'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                      (nsfile, ewfile, udfile, bbpfile, self.log))
        bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

        # Copy output velocity file to output dir
        shutil.copy2(bbpfile, a_outdir)

        # Now, create the acceleration file by combining the 3 acc components

        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "temp_%s.amp.acc.000" % (sta_base))
        ewfile = os.path.join(a_tmpdir,
                              "temp_%s.amp.acc.090" % (sta_base))
        udfile = os.path.join(a_tmpdir,
                              "temp_%s.amp.acc.ver" % (sta_base))
        bbpfile = os.path.join(a_tmpdir,
                               "%d.%s.acc.bbp" % (sim_id, site))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        progstring = ("%s wcc2bbp=1 units=cm/s/s " %
                      (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                      'title="Sim NGAH, stat=%s" ' % (site) +
                      'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                      (nsfile, ewfile, udfile, bbpfile, self.log))
        bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

        # Copy output acceleration file to output dir
        shutil.copy2(bbpfile, a_outdir)

    def run(self):
        """
        Run the GP WccSiteamp 2014 module
        """
        print("GP Site Response".center(80, '-'))

        self.install = InstallCfg.getInstance()
        install = self.install
        self.config = WccSiteampCfg(self.vmodel_name, self.method)
        config = self.config

        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.wcc_siteamp_%s.log" % (sim_id, sta_base))

        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))

        progstring = "mkdir -p %s" % (a_tmpdir)
        bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

        #
        # Read and parse the station list with this call
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for sites in site_list:
            site = sites.scode
            vs30 = sites.vs30
            if vs30 > config.VREF_MAX:
                vs30 = config.VREF_MAX

            print("*** WccSiteamp Processing station %s..." % (site))

            if self.method == "GP":
                self.process_separate_seismograms(site, sta_base,
                                                  vs30, a_tmpdir)
            elif self.method == "SDSU" or self.method == "EXSIM" or self.method == "UCSB":
                self.process_hybrid_seismogram(site, sta_base, vs30,
                                               a_tmpdir, a_outdir)

        print("GP Site Response Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename((sys.argv[0]))))
    ME = WccSiteamp(sys.argv[1], sys.argv[2],
                    sys.argv[3], sim_id=int(sys.argv[4]))
    ME.run()
