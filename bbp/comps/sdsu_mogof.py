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

Broadband Platform Version of SDSU MO-GOF
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil
import numpy as np
from decimal import Decimal, getcontext

# Import Broadband modules
import bband_utils
import plot_seismograms as PS
import validation_cfg
from sdsu_mogof_cfg import SDSUMOGofCfg
from install_cfg import InstallCfg
from station_list import StationList
import plot_value_map
import fault_utils

class SDSUMOGoF(object):
    """
    Implement SDSU Mayhew-Olson GOF as a Broadband Component
    """

    def __init__(self, i_r_stations, i_weights, plot_map,
                 i_a_datadir, i_format,
                 i_comparison_label, sim_id=0):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.gof_weights = i_weights
        self.a_datadir = i_a_datadir
        self.format = i_format
        self.comp_label = i_comparison_label

        # Compute GOF, one station at a time
        self.single_stat_run = True
        self.num_stations = 0
        self.num_timesteps = 0
        self.len_seismo = 0.0
        self.input_set_1 = []
        self.input_set_2 = []
        self.config = SDSUMOGofCfg()
        self.config.cfggof["weights"] = self.gof_weights
        self.config.cfggof["input_param"] = self.format
        self.install = InstallCfg.getInstance()
        self.log = os.path.join(self.install.A_OUT_LOG_DIR,
                                str(self.sim_id),
                                "%d.sdsu_mogof.log" % (self.sim_id))
        self.site_fit = []
        self.metric_count = 0
        self.plot_gof_map = plot_map
        self.hypo = []

    def open_file(self, file_name, mode):
        try:
            file = open(file_name, mode)
        except IOError as err:
            print("ERROR (sdsu_mogof: Unable to open the file",
                  file_name, "Ending program.\n", err)
            sys.exit(-1)
        else:
            return file

    def get_file_list(self, filename):
        files = None

        try:
            fp = open(filename, 'r')
            files = fp.readlines()
            for i in range(0, len(files)):
                files[i] = os.path.abspath(files[i].strip())
            fp.close()
        except:
            print("ERROR (sdsu_mogof): Failed to read in files from %s" %
                  (filename))
            sys.exit(-1)
        else:
            return files

    def check_input_list(self, input_set):

        for seismo in input_set:

            print("Checking seismogram file %s" % (seismo))
            if not os.path.exists(seismo):
                print("ERROR (sdsu_mogof): Seismogram %s does not exist" %
                      (seismo))
                sys.exit(-1)

            sp = open(seismo, 'r')
            samples = sp.readlines()
            sp.close()
            header = 0

            # Check that the BBP file is formatted correctly
            num_samples = len(samples)
            tstart = -1.0
            while samples[num_samples - 1] == "":
                num_samples = num_samples - 1
            for i in range(0, num_samples):
                tokens = samples[i].split()
                if tokens[0] != '#':
                    if len(tokens) != 4:
                        print("ERROR (sdsu_mogof): Seismogram "
                              "%s is incorrectly formatted" % (seismo))
                        print(tokens)
                        sys.exit(-1)
                    if tstart == -1.0:
                        tstart = float(tokens[0])
                else:
                    header += 1

            num_samples -= header
            # Check number of timesteps
            if self.num_timesteps == 0 or self.single_stat_run:
                self.num_timesteps = num_samples
                print("Determined num_timesteps = %d" % (self.num_timesteps))
            else:
                if num_samples != self.num_timesteps:
                    print("ERROR (sdsu_mogof): Found %d samples in %s, expected %d" %
                          (num_samples, seismo, self.num_timesteps))
                    sys.exit(-1)

            # Save length of seismogram
            if self.len_seismo < 1.0 or self.single_stat_run:
                self.len_seismo = float(samples[num_samples - 1].split()[0]) - tstart
                print("Determined len_seismo = %f s" % (self.len_seismo))

        return 0

    def concat_input_set(self, input_set, concat_filename):

        try:
            # Open concatenated file
            print("Creating input file set %s" % (concat_filename))
            cp = open(concat_filename, 'w')

            for input in input_set:
                print("Concatenating input file %s" % (input))
                # Read in file
                input_file = open(input, 'r')
                lines = input_file.readlines()
                input_file.close()

                # Append to concatenated file
                for line in lines:
                    # Ensure newline is present
                    if line.strip()[0] != '#' and line.strip()[0] != "\n":
                        if line[-1:] != "\n":
                            line = line + "\n"
                        cp.write(line)
            cp.close()
        except:
            print("ERROR (sdsu_mogof): Failed to concatenate input set")
            sys.exit(-1)

        return 0

    def write_param_dat(self, paramfile):

        try:
            pp = open(paramfile, 'w')
            pp.write("%s\n" % (self.config.INPUT_SEISMO_1))
            pp.write("%s\n" % (self.config.INPUT_SEISMO_2))
            pp.write("%d\n" % (0))
#            pp.write("%d\n" % (self.config.cfggof["num_station"]))
            pp.write("%d\n" % (self.num_stations))
            pp.write("%d\n" % (self.num_timesteps))
            pp.write("%d\n" % (self.num_timesteps))
            pp.write("%s\n" % (self.config.cfggof["input_param"]))
            pp.write("%f\n" % (self.len_seismo))

            if self.config.cfggof["use_nga"] == "true":
                pp.write("%s\n" % (self.config.cfgnga["source_mag"]))
                pp.write("%s\n" % (self.config.cfgnga["dip"]))
                pp.write("%s\n" % (self.config.cfgnga["rake"]))
                pp.write("%s\n" % (self.config.cfgnga["depth_coseismic"]))
                pp.write("%s\n" % \
                             (os.path.split(self.config.cfgnga["site_file"])[1]))

            pp.write("%s\n" % (self.config.cfggof["low_cut"]))
            pp.write("%s\n" % (self.config.cfggof["high_cut"]))
#            for n in range(0, 12):
#                pp.write("%f\n" % (1.0))
            pp.write("%s\n" % (self.config.cfggof["weights"]["pga"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["pgv"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["pgd"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["psa"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["spectral_Fit"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["cumulative_energy_fit"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["inelastic_elastic_fit"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["sepctral_acc"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["spec_duration"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["data_energy_release_duration"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["cross_correlation"]))
            pp.write("%s\n" % (self.config.cfggof["weights"]["fourier_spectrum"]))
            pp.close()
        except IOError:
            print("ERROR (sdsu_mogof): Failed to write %s" % (paramfile))
            sys.exit(-1)

        return 0

#    def BBPReader(self, bbpfile):
#        cp = self.open_file(bbpfile, 'rb')
#        lines = cp.readlines()
#        cp.close()
#        dt = 0.0
#        iindex = 1
#        dcount = 0
#        hdr = 0
#        samples = []
#
#        for line in lines:
#            tokens = line.strip().split()
#            if (tokens[0] != '#'):
#                if (len(tokens) != 4):
#                    print "ERROR (sdsu_mogof): Seismogram %s is incorrectly formatted" % (bbpfile)
#                    print tokens
#                    sys.exit(-1)
#                else:
#                    sample = (float(tokens[0]), float(tokens[1]),
#                              float(tokens[2]), float(tokens[3]))
#                    samples.append(sample)
#                    dcount += 1
#            else:
#                hdr += 1
#
#        if dcount > 2:
#            dt = samples[1][0] - timesteps[0][0]
#
#        return samples, dt, dcount

    def bbp_writer(self, bbpfile, samples):
        if samples is None:
            print("ERROR (sdsu_mogof): Need to specify sample list")
            sys.exit(-1)

        # Open file
        bbp = open(bbpfile, 'w')
        for sample in samples:
            bbp.write("%10.6f %16.8G %16.8G %16.8G\n" %
                      (float(sample[0]), float(sample[1]),
                       float(sample[2]), float(sample[3])))
        # Close file
        bbp.close()

        # All done
        return 0

    def get_sample_data(self, bbpfile):
        """
        This function reads a BBP file and determines its DT, it also
        returns the number of samples in dcount
        """
        dt = Decimal(0)
        dcount = 0
        bbp_file = self.open_file(bbpfile, 'r')
        for line in bbp_file:
            line = line.strip()
            if line.startswith("#") or line.startswith("%"):
                # Skip comments
                continue
            tokens = line.split()
            if len(tokens) != 4:
                print("ERROR (sdsu_mogof): "
                      "Seismogram %s is incorrectly formatted" % (bbpfile))
                print(tokens)
                sys.exit(-1)

            dcount += 1
            if dcount < 2:
                dt = Decimal(tokens[0])
            elif dcount == 2:
                dt = Decimal(tokens[0]) - dt

        # Close file
        bbp_file.close()

        # All done!
        return dt, dcount

    def match_sample_rate(self, bbpfile, newdt):
        work_dir = os.getcwd()
        bbpfile_name = os.path.join(work_dir,
                                    os.path.splitext(os.path.basename(bbpfile))[0])
        cmd = ("%s/wcc2bbp " % (self.install.A_GP_BIN_DIR) +
               "nsfile=%s.000 ewfile=%s.090 udfile=%s.ver " %
               (bbpfile_name, bbpfile_name, bbpfile_name) +
               "wcc2bbp=0 < %s >> %s 2>&1" % (bbpfile, self.log))
        bband_utils.runprog(cmd)

        for comp in self.config.COMPS:
            infile = bbpfile_name + "." + comp
            outfile = bbpfile_name + ".ms." + comp
            cmd = ("%s/wcc_resamp_arbdt " % (self.install.A_GP_BIN_DIR) +
                   "newdt=%f infile=%s outfile=%s >> %s 2>&1" %
                   (newdt, infile, outfile, self.log))
            bband_utils.runprog(cmd)

        cmd = ("%s/wcc2bbp wcc2bbp=1 " % self.install.A_GP_BIN_DIR +
               'title="MOGOF" '
               'nsfile=%s.ms.000 ewfile=%s.ms.090 udfile=%s.ms.ver > %s.ms.bbp 2>> %s' %
               (bbpfile_name, bbpfile_name, bbpfile_name, bbpfile_name, self.log))
        bband_utils.runprog(cmd)

        os.remove("%s.000" % bbpfile_name)
        os.remove("%s.090" % bbpfile_name)
        os.remove("%s.ver" % bbpfile_name)
        os.remove("%s.ms.000" % bbpfile_name)
        os.remove("%s.ms.090" % bbpfile_name)
        os.remove("%s.ms.ver" % bbpfile_name)


        return "%s.ms.bbp" % bbpfile_name

    def MatchSampleLength(self, bbpfile, length, match_method="pad"):
        work_dir = os.getcwd()
        fname = os.path.join(work_dir,
                             "%s.ml.bbp" %
                             (os.path.splitext(os.path.basename(bbpfile))[0]))

        shutil.copy2(bbpfile, fname)
        bbpfile = fname
        bfile = self.open_file(bbpfile, 'rb')
        lines = bfile.readlines()
        bfile.close()

        if match_method == "pad":
            tokens = lines[len(lines)-1].strip().split()
            bfile = self.open_file(bbpfile, 'a')
            dt = 0.0
            timestep = 0.0

            if len(tokens) != 4:
                print("ERROR (sdsu_mogof): " +
                      "Seismogram %s is incorrectly formatted" % (bbpfile))
                print(tokens)
                sys.exit(-1)
            else:
                timestep = float(tokens[0])
                dt = timestep - float(lines[len(lines) - 2].strip().split()[0])
                for _ in range(0, length):
                    timestep = timestep + dt
                    bfile.write("%10.6f %16.8G %16.8G %16.8G\n" %
                                (timestep, float(tokens[1]),
                                 float(tokens[2]), float(tokens[3])))
                bfile.close()
        elif match_method == "trim":
            print("ERROR (sdsu_mogof): Seismogram trim feature is unavailable")
            sys.exit(-1)
        return bbpfile

    def match_seis(self, syn_bbp, obs_bbp, length=0):
        syn_dt = Decimal(0)
        syn_dcount = 0
        obs_dt = Decimal(0)
        obs_dcount = 0

        if os.path.exists(syn_bbp):
            syn_dt, syn_dcount = self.get_sample_data(syn_bbp)
        else:
            print("ERROR (sdsu_mogof): Seismogram file %s not found!" %
                  (syn_bbp))
            sys.exit(-1)

        if os.path.exists(obs_bbp):
            obs_dt, obs_dcount = self.get_sample_data(obs_bbp)
        else:
            print("ERROR (sdsu_mogof): Seismogram file %s not found!" %
                  (obs_bbp))
            sys.exit(-1)

        print("Matching Seismograms Syn: "
              "%s (dt: %f, SL: %d) with Obs: %s (dt: %f, SL: %d)" %
              (syn_bbp, syn_dt, syn_dcount, obs_bbp, obs_dt, obs_dcount))
        # Match Sample Rates if required
        # Always try to Upsample
        if syn_dt != obs_dt:
            if syn_dt != Decimal(0) and obs_dt != Decimal(0):
                if obs_dt > syn_dt:
                    # print "SDSU_MOGof: Upsampling Observed Seismogram %s from %f to %f" % (obs_bbp, float(obs_dt), float(syn_dt))
                    obs_bbp = self.match_sample_rate(obs_bbp, float(syn_dt))
                    obs_dt, obs_dcount = self.get_sample_data(obs_bbp)
                else:
                    # print "SDSU_MOGof: Upsampling Synthetic Seismogram %s from %f to %f" % (obs_bbp, float(syn_dt), float(obs_dt))
                    syn_bbp = self.match_sample_rate(syn_bbp, float(obs_dt))
                    syn_dt, syn_dcount = self.get_sample_data(syn_bbp)
            else:
                print("ERROR (sdsu_mogof): Time Step is 0! "
                      "Obs-dt: %f, Syn-dt: %f" % (float(obs_dt), float(syn_dt)))
                sys.exit(-1)

        # Match Record Lengths
        if syn_dcount != obs_dcount or length != 0:
            if length == 0:
                # Lenght was not specified, pad the shorter record
                # with zero values at the end of the record
                length = obs_dcount - syn_dcount
                if length > 0:
                    # print "Matching synthetic seismogram length to %d" % obs_dcount
                    syn_bbp = self.MatchSampleLength(syn_bbp, length, "pad")
                elif length < 0:
                    # print "Matching observed seismogram length to %d" % syn_dcount
                    obs_bbp = self.MatchSampleLength(obs_bbp,
                                                     abs(length), "pad")
            elif length > 0:
                # Length was specified adjust the length of the
                # seismogram to specified length
                if length > syn_dcount:
                    syn_bbp = self.MatchSampleLength(syn_bbp, length, "pad")
                elif length < syn_dcount:
                    syn_bbp = self.MatchSampleLength(syn_bbp, length, "trim")
                if length > obs_dcount:
                    obs_bbp = self.MatchSampleLength(obs_bbp, length, "pad")
                elif length < obs_dcount:
                    obs_bbp = self.MatchSampleLength(obs_bbp, length, "trim")

        return syn_bbp, obs_bbp

    def summarize_results(self, gof_results, summ_results):
        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(self.sim_id))
        if self.config.cfggof["num_station"] >= 3 and self.plot_gof_map:
            # Check if a station file with subset of stations exists and use that
            stat_file = "%d_stats.txt" % (self.sim_id)
            if os.path.isfile(os.path.join(a_indir, stat_file)):
                pvm = plot_value_map.PlotValueMap(stat_file,
                                                  self.sim_id, self.hypo)
            else:
                pvm = plot_value_map.PlotValueMap(self.r_stations,
                                                  self.sim_id, self.hypo)
        else:
            if self.plot_gof_map:
                # Disable generating the GOF map with less than 3 stations
                print("Note (SDSU - GOF): Less than 3 stations found in " +
                      "station list (%d), GOF map plots will be disabled!" %
                      self.config.cfggof["num_station"])
                self.plot_gof_map = False

        self.site_fit = np.array([[0.0 for i in range(4)] for j in range(self.config.cfggof["num_station"])])
        self.metric_count = 0
        # Compute site_fit
        if self.config.cfggof["weights"]["pga"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["pga"])
            PGA_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (PGA_val * float(self.config.cfggof["weights"]["pga"])))
            self.metric_count += 1
            if self.format == 'A':
                self.plotOverlays(self.input_set_1, self.input_set_2, PGA_val)
            if self.plot_gof_map:
                pvm.run(PGA_val, "PGA")

        if self.config.cfggof["weights"]["pgv"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["pgv"])
            PGV_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (PGV_val * self.config.cfggof["weights"]["pgv"]))
            self.metric_count += 1
            if self.format == 'V':
                self.plotOverlays(self.input_set_1, self.input_set_2, PGV_val)
            if self.plot_gof_map:
                pvm.run(PGV_val, "PGV")

        if self.config.cfggof["weights"]["pgd"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["pgd"])
            PGD_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["pgd"] * PGD_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(PGD_val, "PGD")

        if self.config.cfggof["weights"]["psa"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["psa"])
            PSA_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["psa"] * PSA_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(PSA_val, "PSA")

        if self.config.cfggof["weights"]["spectral_Fit"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["spectral_Fit"])
            SpFit_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["spectral_Fit"]
                              * SpFit_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(SpFit_val, "Spectral Fit")

        if self.config.cfggof["weights"]["cumulative_energy_fit"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["cumulative_energy_fit"])
            EnFit_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["cumulative_energy_fit"]
                              * EnFit_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(EnFit_val, "Cumulative Energy Fit")

        if self.config.cfggof["weights"]["data_energy_release_duration"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["data_energy_release_duration"])
            Dur_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["data_energy_release_duration"]
                              * Dur_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(Dur_val, "Energy release Duration")

        if self.config.cfggof["weights"]["inelastic_elastic_fit"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["inelastic_elastic_fit"])
            InElEl_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["inelastic_elastic_fit"]
                              * InElEl_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(InElEl_val, "In-elastic - Elastic Fit")

        if self.config.cfggof["weights"]["sepctral_acc"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["sepctral_acc"])
            SAFit_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["sepctral_acc"]
                              * SAFit_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(SAFit_val, "Spectral Acceleration")

        if self.config.cfggof["weights"]["spec_duration"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["spec_duration"])
            SpecDur_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["spec_duration"]
                              * SpecDur_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(SpecDur_val, "Spectral Duration")

        if self.config.cfggof["weights"]["cross_correlation"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["cross_correlation"])
            CCFit_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["cross_correlation"]
                              * CCFit_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(CCFit_val, "Cross Correlation")

        if self.config.cfggof["weights"]["fourier_spectrum"] > 0.0:
            fname = os.path.join(self.config.cfggof["output_dir"],
                                 self.config.cfggof["file"]["fourier_spectrum"])
            FSComp_val = np.array(self.get_gof_metric(fname, start_col=0))
            self.site_fit = (self.site_fit +
                             (self.config.cfggof["weights"]["fourier_spectrum"]
                              * FSComp_val))
            self.metric_count += 1
            if self.plot_gof_map:
                pvm.run(FSComp_val, "Fourier Spectrum")

        #Calculate site_fit Values
        self.site_fit = self.site_fit / float(self.metric_count)
        gof_results = os.path.join(self.config.cfggof["output_dir"],
                                   "GOF.list")
        np.savetxt(gof_results, self.site_fit, delimiter='\t')
        if self.plot_gof_map:
            pvm.run(self.site_fit, "Site Fit")

        # Read in the file
        if not os.path.exists(gof_results):
            print("ERROR (sdsu_mogof): GOF Results file was not found!")
            sys.exit(-1)
        gp = self.open_file(gof_results, 'r')
        lines = gp.readlines()
        gp.close()

        fname = os.path.join(self.config.cfggof["output_dir"],
                             "stat_comp.txt")
        if not os.path.exists(fname):
            print("ERROR (sdsu_mogof): GOF Stations file was not found!")
            sys.exit(-1)
        comp_stat_file = self.open_file(fname, 'r')
        comp_stats = comp_stat_file.readlines()
        comp_stat_file.close()

        sp = open(summ_results, 'w')
        header = "%-24s\t%-24s\tOverall\tX\tY\tZ" % ("Sta (input 1)", "Sta (input 2)")
        sep = "=" * 96
        print("\nResults:")
        print(header)
        print(sep)
        sp.write("#%s\n" % (header))
        sp.write("#%s\n" % (sep))
        stnum = 0
        for line in lines:
            stripped = line.strip()
            if stripped[0] != '%':
                tokens = stripped.split()
                if stnum >= self.config.cfggof["num_station"]:
                    print("More results reported than input stations!")
                    return 1

                stats = comp_stats[stnum].strip().split()
                input_1 = stats[0]
                input_2 = stats[1]
                overall = float(tokens[0])
                x = float(tokens[1])
                y = float(tokens[2])
                z = float(tokens[3])
                results = ("%-24s\t%-24s\t%6.2f\t%6.2f\t%6.2f\t%6.2f" %
                           (input_1, input_2, overall, x, y, z))
                print(results)
                sp.write("%s\n" % (results))
                stnum += 1

        sp.close()

        print("")
        print("Saved summary in %s" % (summ_results))
        return 0

    def get_gof_metric(self, filename, start_col=0):
        gof_metric = []
        metric_file = self.open_file(filename, 'r')
        m_data = metric_file.readlines()

        for line in m_data:
            if line.strip().startswith('#') or line.strip().startswith('%'):
                continue
            else:
                pieces = line.split()
                if start_col == 0:
                    if len(pieces) == 7:
                        # slon, slat, site, overall fit, x fit, y fit, z fit
                        start_col = 3
                gof_metric.append([float(s) for s in pieces[start_col:]])
        return gof_metric

    #def set_gof_metric(self, filename, m_data, stat=None, start_col=0):
    def set_gof_metric(self, filename, m_data, stat=None):
        metric_file = self.open_file(filename, 'a')
        sample = ""
        for line in m_data:
            sample += "%-10s" % "\t".join(map(str, line))
            if stat is not None:
                sample = "%-8s\t%-8s\t%-12s\t%s" % (stat.lon, stat.lat,
                                                    stat.scode, sample)
            metric_file.write("%s\n" % sample)
        metric_file.close()
        return

    def plotOverlays(self, filelist_syn, filelist_obs, gofdata):

        if len(filelist_syn) != len(filelist_obs):
            print("ERROR (sdsu_mogof): "
                  "Number of stations in sets 1 and 2 differ!")
            sys.exit(-1)
        if len(filelist_syn) != len(gofdata):
            print("ERROR (sdsu_mogof): Number of stations "
                  "%d is not equal to lenght of GOF data %d!" %
                  (len(filelist_syn), len(gofdata)))
            sys.exit(-1)

        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR,
                                str(self.sim_id),
                                "gof_plots")
        if not os.path.exists(a_outdir):
            os.mkdir(a_outdir)

            if self.format == "V":
                ylabel = "Velocity (cm/s)"
                goflabel = ["GOF", "PGV"]
            elif self.format == "A":
                ylabel = "Acceleration (cm/s/s)"
                goflabel = ["GOF", "PGA"]

            for i in range(0, len(filelist_syn)):
                obs_bbp = filelist_obs[i]
                syn_bbp = filelist_syn[i]
                gdata = gofdata[i]
                stat = os.path.basename(obs_bbp)
                stat = stat.split(".")[0]

                outfile = os.path.join(a_outdir, "%s_%s_match_%s_overlay.png" %
                                       (self.sim_id, stat, self.format))
                print("PS.plot_overlay(%s, %s, %s, Observed, Run %s," %
                      (stat, obs_bbp, syn_bbp, self.sim_id) +
                      "%s, %s, %r, %r" %
                      (outfile, ylabel, goflabel, gdata))
                PS.plot_overlay(stat, obs_bbp, syn_bbp, "Observed",
                                "Run %s" % self.sim_id, outfile,
                                ylabel, goflabel, gdata)

            return

    def do_gof(self, filelist_syn, filelist_obs, stat=None):

        sim_id = self.sim_id
        run_dir = self.config.cfggof["work_dir"]

        self.config.cfggof["input_set_1"] = os.path.join(run_dir,
                                                         "%d_syn_list.input" %
                                                         (sim_id))
        self.config.cfggof["input_set_2"] = os.path.join(run_dir,
                                                         "%d_obs_list.input" %
                                                         (sim_id))

        inpsfile = self.open_file(self.config.cfggof["input_set_1"], 'w')
#        num_stations = 0
        for filename in filelist_syn:
            inpsfile.write("%s\n" % (filename))
#            num_stations += 1
        inpsfile.close()
        inpofile = self.open_file(self.config.cfggof["input_set_2"], 'w')
        for filename in filelist_obs:
            inpofile.write("%s\n" % (filename))
        inpofile.close()

        #Concatenate input BBP files into single GOF
        # User data in a_outdir, reference data in self.install.A_GF_DIR/gp/data/AccBBP
        #  - Check input files exist
        if not os.path.exists(self.config.cfggof["input_set_1"]):
            print("ERROR (sdsu_mogof): Input_set_1 file %s does not exist" %
                  (self.config.cfggof["input_set_1"]))
            sys.exit(-1)
        if not os.path.exists(self.config.cfggof["input_set_2"]):
            print("ERROR (sdsu_mogof): Input_set_2 file %s does not exist" %
                  (self.config.cfggof["input_set_2"]))
            sys.exit(-1)

        # Get individual input file lists
        input_set_1 = self.get_file_list(self.config.cfggof["input_set_1"])
        input_set_2 = self.get_file_list(self.config.cfggof["input_set_2"])
        if input_set_1 is None or input_set_2 is None:
            print("ERROR (sdsu_mogof): Empty file lists found")
            sys.exit(-1)

        # Check each file in each filelist
        if len(input_set_1) != len(input_set_2):
            print("ERROR (sdsu_mogof): "
                  "Number of stations in sets 1 and 2 differ")
            sys.exit(-1)
        if self.check_input_list(input_set_1) != 0:
            print("ERROR (sdsu_mogof): Input_set_1 is not valid")
            sys.exit(-1)
        if self.check_input_list(input_set_2) != 0:
            print("ERROR (sdsu_mogof): Input_set_2 is not valid")
            sys.exit(-1)

        # Concatenate each input set into a single GOF formatted file
        # and save to work_dir
        concat_file_1 = os.path.join(self.config.cfggof["work_dir"],
                                     self.config.INPUT_SEISMO_1)
        concat_file_2 = os.path.join(self.config.cfggof["work_dir"],
                                     self.config.INPUT_SEISMO_2)
        if self.concat_input_set(input_set_1, concat_file_1) != 0:
            print("ERROR (sdsu_mogof): Failed to concatenate input set 1")
            sys.exit(-1)
        if self.concat_input_set(input_set_2, concat_file_2) != 0:
            print("ERROR (sdsu_mogof): Failed to concatenate input set 2")
            sys.exit(-1)

#        self.config.cfggof["num_station"]= len(input_set_1)
        self.config.cfggof["timesteps_set_1"] = self.num_timesteps
        self.config.cfggof["timesteps_set_2"] = self.num_timesteps
        self.config.cfggof["seismo_length"] = self.len_seismo

        # If NGA, copy site_file to work_dir
        if self.config.cfggof["use_nga"]:
            print("ERROR (sdsu_mogof): "
                  "USE_NGA option has not been implemented. Aborting...")
            #shutil.copy(self.config.cfgnga["site_file"], self.config.cfggof["work_dir"])
            sys.exit(-1)

        # Construct parameter file in working_dir   <<<<<<
        paramfile = os.path.join(self.config.cfggof["work_dir"],
                                 self.config.PARAM_DAT_FILE)
        if self.write_param_dat(paramfile) != 0:
            print("ERROR (sdsu_mogof): Failed to write %s" % (paramfile))
            sys.exit(-1)

        # Save path to GOF output dir
        gof_out_dir = os.path.abspath(os.path.join(self.config.cfggof["work_dir"], "out"))
        # Check if "out" subdir exists
        if os.path.exists(gof_out_dir):
            # Remove "out" subdir
            shutil.rmtree(gof_out_dir)

        # (Re)Create "out" subdir
        os.mkdir(gof_out_dir)

        # Execute GOF codes

        # Run the code
        if self.config.cfggof["use_nga"] == "true":
            # This has not been implemented!!!!
            # Just abort (there's already a check about 30 lines above)
            sys.exit(-1)
            # cmd = ['%s' % (gof_nga_bin)]
            # print "Executing cmd: %s" % (str(cmd))
            # bband_utils.runprog(cmd)
        else:
#            cmd = ['%s' % (gof_bin)]
#            print "Executing cmd: %s" % (str(cmd))
#            bband_utils.runprog(cmd)
#
#            site_fit = np.array([[0.0 for i in range(4)] for j in range(num_stations)])
            self.metric_count = 0

            PGX_run = False

            if self.config.cfggof["weights"]["pga"] > 0.0:
                PGX_run = True
                cmd = '%s >> %s 2>&1' % (self.config.GOF_PGX_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["pga"])
                PGA_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["pga"])
                self.set_gof_metric(fname, PGA_val, stat)
#                site_fit = site_fit + (PGA_val * np.float(self.config.cfggof["weights"]["pga"]))
                self.metric_count += 1

            if self.config.cfggof["weights"]["pgv"] > 0.0:
                if PGX_run == False:
                    PGX_run = True
                    cmd = '%s >> %s 2>&1' % (self.config.GOF_PGX_BIN, self.log)
                    print("Executing cmd: %s" % (str(cmd)))
                    bband_utils.runprog(cmd)

                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["pgv"])
                PGV_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["pgv"])
                self.set_gof_metric(fname, PGV_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["pgv"] * PGV_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["pgd"] > 0.0:
                if PGX_run == False:
                    PGX_run = True
                    cmd = '%s >> %s 2>&1' % (self.config.GOF_PGX_BIN, self.log)
                    print("Executing cmd: %s" % (str(cmd)))
                    bband_utils.runprog(cmd)

                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["pgd"])
                PGD_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["pgd"])
                self.set_gof_metric(fname, PGD_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["pgd"] * PGD_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["psa"] > 0.0:
                if PGX_run == False:
                    PGX_run = True
                    cmd = '%s >> %s 2>&1' % (self.config.GOF_PGX_BIN, self.log)
                    print("Executing cmd: %s" % (str(cmd)))
                    bband_utils.runprog(cmd)

                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["psa"])
                PSA_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["psa"])
                self.set_gof_metric(fname, PSA_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["psa"] * PSA_val)
                self.metric_count += 1


            if self.config.cfggof["weights"]["spectral_Fit"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_SpFit_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["spectral_Fit"])
                SpFit_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["spectral_Fit"])
                self.set_gof_metric(fname, SpFit_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["spectral_Fit"] * SpFit_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["cumulative_energy_fit"] > 0.0 \
               or self.config.cfggof["weights"]["data_energy_release_duration"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_DCumEn_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                if self.config.cfggof["weights"]["cumulative_energy_fit"] > 0.0:
                    fname = os.path.join(self.config.cfggof["work_dir"],
                                         "out",
                                         self.config.cfggof["file"]["cumulative_energy_fit"])
                    EnFit_val = np.array(self.get_gof_metric(fname, start_col=0))
                    fname = os.path.join(self.config.cfggof["output_dir"],
                                         self.config.cfggof["file"]["cumulative_energy_fit"])
                    self.set_gof_metric(fname, EnFit_val, stat)
#                    site_fit = site_fit + (self.config.cfggof["weights"]["cumulative_energy_fit"] * EnFit_val)
                    self.metric_count += 1
                if self.config.cfggof["weights"]["data_energy_release_duration"] > 0.0:
                    fname = os.path.join(self.config.cfggof["work_dir"],
                                         "out",
                                         self.config.cfggof["file"]["data_energy_release_duration"])
                    Dur_val = np.array(self.get_gof_metric(fname, start_col=0))
                    fname = os.path.join(self.config.cfggof["output_dir"],
                                         self.config.cfggof["file"]["data_energy_release_duration"])
                    self.set_gof_metric(fname, Dur_val, stat)
#                    site_fit = site_fit + (self.config.cfggof["weights"]["data_energy_release_duration"] * Dur_val)
                    self.metric_count += 1

            if self.config.cfggof["weights"]["inelastic_elastic_fit"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_InElFit_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["inelastic_elastic_fit"])
                InElEl_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["inelastic_elastic_fit"])
                self.set_gof_metric(fname, InElEl_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["inelastic_elastic_fit"] * InElEl_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["sepctral_acc"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_SAFit16_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["sepctral_acc"])
                SAFit_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["sepctral_acc"])
                self.set_gof_metric(fname, SAFit_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["sepctral_acc"] * SAFit_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["spec_duration"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_SpecDurFit_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["spec_duration"])
                SpecDur_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["spec_duration"])
                self.set_gof_metric(fname, SpecDur_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["spec_duration"] * SpecDur_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["cross_correlation"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_CCFit_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["cross_correlation"])
                CCFit_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["cross_correlation"])
                self.set_gof_metric(fname, CCFit_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["cross_correlation"] * CCFit_val)
                self.metric_count += 1

            if self.config.cfggof["weights"]["fourier_spectrum"] > 0.0:
                cmd = '%s >> %s 2>&1' % (self.config.GOF_FSComp_BIN, self.log)
                print("Executing cmd: %s" % (str(cmd)))
                bband_utils.runprog(cmd)
                fname = os.path.join(self.config.cfggof["work_dir"],
                                     "out",
                                     self.config.cfggof["file"]["fourier_spectrum"])
                FSComp_val = np.array(self.get_gof_metric(fname, start_col=0))
                fname = os.path.join(self.config.cfggof["output_dir"],
                                     self.config.cfggof["file"]["fourier_spectrum"])
                self.set_gof_metric(fname, FSComp_val, stat)
#                site_fit = site_fit + (self.config.cfggof["weights"]["fourier_spectrum"] * FSComp_val)
                self.metric_count += 1

            #Calculate site_fit Values

#            site_fit = site_fit / float(self.metric_count)
#            fname = "%s/out/GOF.list" % (self.config.cfggof["work_dir"])
#            np.savetxt(fname, site_fit, delimiter='\t')
        return

    def run(self):
        """
        Runs the SDSU MoGOF code
        """
        print("SDSU MOGoF".center(80, '-'))

        getcontext().prec = 10
        # Required Inputs: sim_id
        sim_id = self.sim_id
        old_cwd = os.getcwd()
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        a_indir = os.path.join(self.install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(self.install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(self.install.A_TMP_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(sim_id))
        a_statfile = os.path.join(self.install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        if not os.path.exists(a_statfile):
            # We need the station file
            print("ERROR (sdsu_mogof): Cannot find station list file %s" %
                  (a_statfile))
            sys.exit(-1)

        self.config.cfggof["output_dir"] = a_outdir

        # Make sure the output and tmp directories exist
        if not os.path.exists(a_tmpdir):
            os.mkdir(a_tmpdir)

        if not os.path.exists(a_outdir):
            os.mkdir(a_outdir)

        if self.plot_gof_map:
            # Let's figure out the hypocenter location
            val_obj = None
            if self.comp_label != "" and self.comp_label is not None:
                # If validation, get validation object
                val_obj = validation_cfg.VE_EVENTS.get_event_by_name(self.comp_label)
            if val_obj is not None:
                # Look to see if we have a srf file specified
                srffile = val_obj.get_input("GP", "srf")
            else:
                srffile = ""

            # If we have SRF file, get hypocenter from it, otherwise
            # look for an srf file in the indata directory
            if srffile != "":
                self.hypo = fault_utils.get_hypocenter(srffile)
            else:
                entries = os.listdir(a_indir)
                candidate_list = []
                for entry in entries:
                    if entry.startswith("xyz_"):
                        # Skip xyz srf file created by bbtoolbox.py
                        continue
                    if entry.endswith(".srf"):
                        candidate_list.append(entry)
                if len(candidate_list) == 1:
                    srffile = os.path.join(a_indir, candidate_list[0])
                    self.hypo = fault_utils.get_hypocenter(srffile)

        # Create a site1d run dir
        run_dir = os.path.join(a_tmpdir, "MOGof")
        if not os.path.exists(run_dir):
            os.mkdir(run_dir)
            self.config.cfggof["work_dir"] = run_dir
        else:
            # - Clean MOGof run dir
            for root, dirs, files in os.walk(run_dir):
                for f in files:
                    os.unlink(os.path.join(root, f))
                for d in dirs:
                    shutil.rmtree(os.path.join(root, d))

        os.chdir(run_dir)

        # Read and parse the station list with this call
        slo = StationList(a_statfile)
        if slo is None:
            print("ERROR (sdsu_mogof): Cannot open station list %s" %
                  (a_statfile))
            sys.exit(-1)

        site_list = slo.get_station_list()

        print("Opening Station list %s." % (a_statfile))

        #Build GOF formated input file from BPP files.
        # - Check if the output BBP file exists for a station in the station list
        # - Concatenate it to end of the input file.
        filelist_syn = []
        filelist_obs = []
        self.num_stations = 0
        self.config.cfggof["num_station"] = 0
        stats = []
        fname = os.path.join(self.config.cfggof["output_dir"], "stat_comp.txt")
        comp_stat_file = self.open_file(fname, 'a')
        print("Running SDSU MO-GOF codes")
        for sites in site_list:
            slon = float(sites.lon)
            slat = float(sites.lat)
            site = sites.scode

            expected_file = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                         (sim_id, site))
            if not os.path.exists(expected_file):
                # Just skip it
                print("Couldn't find file %s. This is not necessarily an error, as you may have run with a subset of a stations. Goodness of fit will continue with available stations." % (expected_file))
                continue

            obs_file = os.path.join(a_tmpdir_seis, "%s.bbp" % (site))
            if not os.path.exists(obs_file):
                # Just skip it
                print("Couldn't find observed seismogram file %s. Goodness of fit will continue with available stations." % (obs_file))
                continue
            stats.append((slon, slat, site))
            self.config.cfggof["num_station"] += 1

            expected_file, obs_file = self.match_seis(expected_file,
                                                      obs_file, 0)
            # Add file to input file list
            filelist_syn.append(expected_file)
            self.input_set_1.append(expected_file)
            filelist_obs.append(obs_file)
            self.input_set_2.append(obs_file)
            comp_stat_file.write("%s\t\t%s\n" %
                                 (os.path.basename(expected_file),
                                  os.path.basename(obs_file)))
            self.num_stations += 1

            if self.single_stat_run:
                self.do_gof(filelist_syn, filelist_obs, stat=sites)
                self.num_stations = 0
                filelist_syn = []
                filelist_obs = []

        if self.config.cfggof["num_station"] != len(site_list):
            print("Writing a subset of stations with valid BBP files for GOF!")
            fname = os.path.join(a_indir, "%d_stats.txt" % (self.sim_id))
            stat_file = self.open_file(fname, 'a')
            for station in stats:
                stat_file.write("%.3f\t%.3f\t%s\n" % (station[0], station[1],
                                                      station[2]))
            stat_file.close()

#        if self.single_stat_run == False:
#            self.do_gof(filelist_syn, filelist_obs)

        comp_stat_file.close()
        # Summarize results
        gof_results = os.path.join(self.config.cfggof["output_dir"], "GOF.list")
        summ_results = os.path.join(self.config.cfggof["output_dir"],
                                    self.config.SUMMARY_FILE)
        if self.summarize_results(gof_results, summ_results) != 0:
            print("ERROR (sdsu_mogof): Failed to summarize results")
            sys.exit(-1)

        os.chdir(old_cwd)
        print("SDSU MOGoF Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = SDSUMOGoF(sys.argv[1], sys.argv[2], sys.argv[3],
                   sys.argv[4], sys.argv[5], sys.argv[6],
                   sys.argv[7], sim_id=int(sys.argv[8]))
    ME.run()
