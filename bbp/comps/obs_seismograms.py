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

Module to prepare observation data files for processing by the
Broadband Platform
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList
import bbp_formatter
import gmpe_config
from rotd50 import RotD50
from correct_psa import CorrectPSA

SUPPORTED_OBS_FORMATS = ["acc_bbp", "acc_peer", "gmpe"]

class ObsSeismograms(object):

    def __init__(self, i_r_stations,
                 i_a_obsdir, i_obs_format, i_obs_corr,
                 sim_id=0):
        """
        Initialize basic parameters for the ObsSeismograms class
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.a_obsdir = i_a_obsdir
        self.obs_format = i_obs_format
        self.obs_corrections = i_obs_corr
        # Make observed seismograms are in a format we can handle
        if i_obs_format not in SUPPORTED_OBS_FORMATS:
            raise bband_utils.ParameterError("Format %s for " %
                                             (self.obs_format) +
                                             "observed seismograms "
                                             "not supported")

    def run(self):
        """
        This function copies the observed seismograms for the stations
        specified in r_stations to a temporary directory inside the
        tmpdata directory and converts them to the format needed by
        the goodness of fitness module
        """
        print("ObsSeismograms".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.obs_seis.log" %
                                (sim_id))

        # Input, tmp, and output directories
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_seis = os.path.join(a_outdir, "obs_seis_%s" % (sta_base))

        #
        # Make sure the output and tmp directories exist
        #
        dirs = [a_tmpdir, a_tmpdir_seis, a_outdir, a_outdir_seis]
        bband_utils.mkdirs(dirs, print_cmd=False)

        # Station file
        a_statfile = os.path.join(a_indir, self.r_stations)
        # List of observed seismogram files
        filelist = os.listdir(self.a_obsdir)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Inialize the CorrectPSA module
        if self.obs_corrections:
            corr_psa = CorrectPSA(self.r_stations,
                                  "rd50",
                                  os.path.join(a_indir,
                                               self.obs_corrections),
                                  a_tmpdir_seis,
                                  self.sim_id)
        else:
            corr_psa = None

        # Go through each station
        for site in site_list:
            slon = float(site.lon)
            slat = float(site.lat)
            stat = site.scode
            print("==> Processing data for station: %s" % (stat))
            # Since we're using the GP station list, make sure the
            # .rd50 for the station exists.  It might not if we ran the
            # validation with a different station list (like UCSB for
            # Landers)
            expected_rd50_file = os.path.join(a_outdir,
                                              "%d.%s.rd50" %
                                              (sim_id, stat))
            if not os.path.exists(expected_rd50_file):
                # just skip it
                print("Couldn't find file %s. " %
                      (expected_rd50_file) +
                      "This is not necessarily an error, as you may have " +
                      "run with a subset of a stations. Continuing " +
                      "with available stations.")
                continue

            # Ok, we have a calculated rd50 file for this station,
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
                print("==> Converting file: %s" % (r_bbp_file))
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
                # print(r_e_peer_file, r_n_peer_file, r_z_peer_file)
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
            elif self.obs_format == "gmpe":
                # GMPE verification packages don't have actual
                # seismograms, so there's nothing we need to do here!
                r_gmpe_file = "%s-gmpe.ri50" % (stat)
                if r_gmpe_file not in filelist:
                    # No gmpe file for this station
                    continue
                # Copy gmpe file to the tmp obs directory and rename
                # it so that it has a rd50 extension
                a_src_gmpe_file = os.path.join(self.a_obsdir, r_gmpe_file)
                a_dst_rd50_file = os.path.join(a_tmpdir_seis, "%s.rd50" %
                                               (stat))
                shutil.copy2(a_src_gmpe_file, a_dst_rd50_file)
                # Create a copy in outdata that averages all gmpes
                a_avg_rd50_file = os.path.join(a_outdir_seis,
                                               "%s.rd50" % (stat))
                gmpe_config.average_gmpe(stat,
                                         a_src_gmpe_file,
                                         a_avg_rd50_file)
                # All done!
                continue
            else:
                raise bband_utils.ParameterError("Format %s for " %
                                                 (self.obs_format) +
                                                 "observed seismograms "
                                                 "not supported")

            # Run RotD50 on this file
            if corr_psa is not None:
                # First calculate rd50 and psa5 files
                print("===> Calculating RotD50 for station: %s" % (stat))
                RotD50.do_rotd50(a_tmpdir_seis, r_e_peer_file,
                                 r_n_peer_file, r_z_peer_file,
                                 "%s-orig.rd50" % (stat),
                                 self.log)

                # Now we need to correct the RotD50 output using the
                # user-supplied correction factors
                print("===> Correcting PSA for station: %s" % (stat))
                corr_psa.correct_station(stat, "rd50")
            else:
                # Use final names for output files
                print("===> Calculating RotD50 for station: %s" % (stat))
                RotD50.do_rotd50(a_tmpdir_seis, r_e_peer_file,
                                 r_n_peer_file, r_z_peer_file,
                                 "%s.rd50" % (stat),
                                 self.log)
            shutil.copy2(os.path.join(a_tmpdir_seis, "%s.rd50" % (stat)),
                         os.path.join(a_outdir_seis, "%s.rd50" % (stat)))

        print("ObsSeismograms Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    if len(sys.argv) < 6:
        print("Usage: %s " % (os.path.basename(sys.argv[0])) +
              "station_file obs_dir obs_format obs_corr_file sim_id")
        sys.exit(1)
    OBS_SEIS = ObsSeismograms(sys.argv[1], sys.argv[2],
                              sys.argv[3], sys.argv[4],
                              int(sys.argv[5]))
    OBS_SEIS.run()
