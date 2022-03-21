#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
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

Generates FAS GoF plot
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import numpy as np

# Import Broadband modules
import bband_utils
import plot_config
from fas_gof_cfg import FASGofCfg
from install_cfg import InstallCfg
from station_list import StationList
from PlotGOF import PlotGoF

# Import Pynga and its utilities
import pynga.utils as putils

T95 = [ 6.3138, 2.9200, 2.3534, 2.1318, 2.0150, 1.9432, 1.8946,
        1.8595, 1.8331, 1.8125, 1.7959, 1.7823, 1.7709, 1.7613,
        1.7531, 1.7459, 1.7396, 1.7341, 1.7291, 1.7247, 1.7207,
        1.7171, 1.7139, 1.7109, 1.7081, 1.7056, 1.7033, 1.7011,
        1.6991, 1.6973, 1.6955, 1.6939, 1.6924, 1.6909, 1.6896,
        1.6883, 1.6871, 1.6860, 1.6849, 1.6839, 1.6829, 1.6820,
        1.6811, 1.6802, 1.6794, 1.6787, 1.6779, 1.6772, 1.6766,
        1.6759, 1.6753, 1.6747, 1.6741, 1.6736, 1.6730]

def read_bbp_dt(bbp_file):
    """
    Reads BBP file and returns the DT
    """
    # Pick DT from these files
    bbp_dt = None
    input_file = open(bbp_file)
    for line in input_file:
        line = line.strip()
        if line.startswith("#") or line.startswith("%"):
            continue
        # Got to first timestamp. Now, pick two consecutive
        # timestamps values
        bbp_t1 = float(line.strip().split()[0])
        bbp_t2 = float(next(input_file).strip().split()[0])
        # Subtract the two times
        bbp_dt = bbp_t2 - bbp_t1
        # All done!
        break
    input_file.close()

    if bbp_dt is None:
        raise bband_utils.ParameterError("Cannot find DT in %s!" %
                                         (bbp_file))

    return bbp_dt

def rewrite_fas_eas_file(fas_input_file, fas_output_file):
    """
    Reads the fas_input_file, and writes its
    content back without the eas column so that
    it can be used by the GoF tools
    """
    input_file = open(fas_input_file, 'r')
    output_file = open(fas_output_file, 'w')
    output_file.write("# Freq(Hz)\t FAS H1 (cm/s)\t FAS H2 (cm/s)\t "
                   "Smoothed EAS (cm/s)\n")
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        if len(pieces) != 5:
            continue
        pieces = [float(piece) for piece in pieces]
        output_file.write("%2.7E\t%2.7E\t%2.7E\t%2.7E\n" %
                          (pieces[0], pieces[1], pieces[2], pieces[4]))

    input_file.close()
    output_file.close()

class FASGof(object):
    """
    This class generates GOF plots for the FAS data
    """

    def __init__(self, i_r_srcfile, i_r_stations,
                 i_comparison_label, method,
                 cutoff=None, sim_id=0):
        """
        Initialize class instance variables
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_stations = i_r_stations
        self.method = method.lower()
        self.comp_label = i_comparison_label
        self.max_cutoff = cutoff
        self.install = None
        self.config = None
        self.src_keys = None

    def resid2uncer_py(self, residfile, fileroot, comp,
                       nstat, min_cdst, max_cdst):
        """
        Python version of GP resid2uncer code
        """
        min_vs30 = -1e+15
        max_vs30 =  1e+15
        min_xcos = -1e+15
        max_xcos =  1e+15
        min_ycos = -1e+15
        max_ycos =  1e+15
        
        input_file = open(residfile, 'r')
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("%"):
                continue
            break
        # Get list of periods
        per = line.split()[13:]
        per = [float(period) for period in per]
        nper = len(per)
        # print("periods in file = %d\n" % (nper))

        nstat_read = 0
        tmin = []
        tmax = []
        resid = []
        nval = np.zeros(len(per), dtype=int)
        bias = np.zeros(len(per))
        sigma = np.zeros(len(per))
        sigma0 = np.zeros(len(per))
        cl90m = np.zeros(len(per))
        cl90p = np.zeros(len(per))
        
        # Continue reading input file
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            pieces = line.split()
            vs30 = float(pieces[6])
            cdst = float(pieces[7])
            xcos = float(pieces[8])
            ycos = float(pieces[9])
            rdcomp = pieces[12]

            if (vs30 >= min_vs30 and
                vs30 <= max_vs30 and
                cdst >= min_cdst and
                cdst <= max_cdst and
                xcos >= min_xcos and
                xcos <= max_xcos and
                ycos >= min_ycos and
                ycos <= max_ycos and
                rdcomp == comp):

                # Read this station
                tmin.append(float(pieces[10]))
                tmax.append(float(pieces[11]))

                pieces = pieces[13:]
                pieces = [float(piece) for piece in pieces]
                resid.append(pieces)
                
                nstat_read = nstat_read + 1
                if nstat_read > nstat:
                    pass

        input_file.close()

        # Completed reading input file
        # print("nstat_read = %d " % (nstat_read))
        self.uncer(nstat_read, nval, tmin, tmax, nper, per,
                   resid, bias, sigma, sigma0, cl90m, cl90p)

        # Write output to files
        bias_file = open("%s.bias" % (fileroot), 'w')
        for c_per, c_bias in zip(per, bias):
            bias_file.write("%13.5e %13.5e\n" % (c_per, c_bias))
        bias_file.close()

        sigma_file = open("%s.sigma" % (fileroot), 'w')
        for c_per, c_sigma in zip(per, sigma):
            sigma_file.write("%13.5e %13.5e\n" % (c_per, c_sigma))
        sigma_file.close()

        sigma0_file = open("%s.sigma0" % (fileroot), 'w')
        for c_per, c_sigma0 in zip(per, sigma0):
            sigma0_file.write("%13.5e %13.5e\n" % (c_per, c_sigma0))
        sigma0_file.close()

        m90_file = open("%s.m90" % (fileroot), 'w')
        for c_per, c_m90 in zip(per, cl90m):
            m90_file.write("%13.5e %13.5e\n" % (c_per, c_m90))
        m90_file.close()

        p90_file = open("%s.p90" % (fileroot), 'w')
        for c_per, c_p90 in zip(per, cl90p):
            p90_file.write("%13.5e %13.5e\n" % (c_per, c_p90))
        p90_file.close()

    def uncer(self, nstat_read, nval, tmin, tmax, nper, per,
              resid, bias, sigma, sigma0, cl90m, cl90p):
        """
        Calculate GoF parameters
        """
        for i in range(0, nper):
            bias[i] = 0.0
            cl90m[i] = 0.0
            nval[i] = 0

        for i in range(0, nstat_read):
            for j in range(0, nper):
                if per[j] >= tmin[i] and per[j] <= tmax[i]:
                    bias[j] = bias[j] + resid[i][j]
                    cl90m[j] = cl90m[j] + resid[i][j] * resid[i][j]
                    nval[j] = nval[j] + 1

        for j in range(0, nper):
            if nval[j] > 1:
                invn = 1.0 / float(nval[j])
                invn1 = 1.0 / float(nval[j] - 1)

                if nval[j] > 56:
                    ttfac = 1.64 * np.sqrt(invn)
                else:
                    ttfac = T95[nval[j]-2] * np.sqrt(invn)

                bias[j] = bias[j] * invn
                sigma[j] = np.sqrt(invn1 * (cl90m[j] - nval[j] * bias[j] * bias[j]))
                sigma0[j] = np.sqrt(invn * (cl90m[j]))
                cl90m[j] = bias[j] - sigma[j] * ttfac
                cl90p[j] = bias[j] + sigma[j] * ttfac
            else:
                bias[j] = 0.0
                sigma[j] = 0.0
                sigma0[j] = 0.0
                cl90m[j] = 0.0
                cl90p[j] = 0.0

    def summarize_fas(self, site_list, a_outdir):
        """
        Summarizes all FAS data and creates the FAS GOF plot
        """
        sim_id = self.sim_id
        install = self.install
        config = self.config

        freq_ranges = plot_config.FAS_GOF_FREQ
        lfreq=freq_ranges[self.method]['freq_low']
        hfreq=freq_ranges[self.method]['freq_high']

        fas_residfile = os.path.join(a_outdir, "%s-%d.fas-resid.txt" %
                                      (self.comp_label, sim_id))
        for comp in config.COMPS_FAS:
            # Build paths and check lengths
            fileroot = os.path.join(a_outdir, "%s-%d_r%d-%d-fas-%s" %
                                    (self.comp_label, sim_id, config.MIN_CDST,
                                     self.max_cutoff, comp))
            bband_utils.check_path_lengths([fas_residfile, fileroot],
                                           bband_utils.GP_MAX_FILENAME)
            
            self.resid2uncer_py(fas_residfile, fileroot, comp,
                                len(site_list), config.MIN_CDST,
                                self.max_cutoff)

            #resid2uncer_bin = os.path.join(install.A_GP_BIN_DIR,
            #                               "resid2uncer_varN")
            #cmd = ("%s " % (resid2uncer_bin) +
            #       "residfile=%s fileroot=%s " % (fas_residfile, fileroot) +
            #       "comp=%s nstat=%d nper=256 " % (comp, len(site_list)) +
            #       "min_cdst=%d max_cdst=%d >> %s 2>&1" %
            #       (config.MIN_CDST, self.max_cutoff, self.log))
            #bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # Plot GOF for FAS data
        plot_mode = 'rd50'
        fileroot = '%s-%d_r0-%d-fas' % (self.comp_label, sim_id, self.max_cutoff)
        plottitle = ("GOF Comparison between %s and simulation %d" %
                     (self.comp_label, sim_id))
        plotter = PlotGoF()
        plotter.plot_fas_gof(plottitle, fileroot, a_outdir,
                             a_outdir, cutoff=self.max_cutoff,
                             lfreq=lfreq, hfreq=hfreq)

    def run(self):
        """
        This function in the main entry point for this module.
        It runs the fas gof component.
        """
        print("FAS GoF".center(80, '-'))

        # Initialize basic variables
        self.install = InstallCfg.getInstance()
        self.config = FASGofCfg()
        install = self.install
        config = self.config
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.fas_gof.log" %
                                (sim_id))

        # Input, tmp, and output directories
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir_fas = os.path.join(a_outdir, "FAS")

        # Source file, parse it!
        a_srcfile = os.path.join(install.A_IN_DATA_DIR,
                                 str(sim_id),
                                 self.r_srcfile)
        self.src_keys = bband_utils.parse_src_file(a_srcfile)

        # Station file
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # check cutoff value
        if self.max_cutoff is None:
            self.max_cutoff = config.MAX_CDST

        print_header_fas = 1
        # Remove fas resid file
        fas_resid_output = os.path.join(a_outdir, "%s-%d.fas-resid.txt" %
                                        (self.comp_label, sim_id))
        if os.path.exists(fas_resid_output):
            os.remove(fas_resid_output)

        for site in site_list:
            slon = float(site.lon)
            slat = float(site.lat)
            stat = site.scode

            # Pick up DT from simulated file
            acc_file = "%d.%s.acc.bbp" % (sim_id, stat)
            input_syn_acc_file = os.path.join(a_outdir, acc_file)
            syn_dt = read_bbp_dt(input_syn_acc_file)
            max_syn_freq = 1.0 / (syn_dt * 2)
            if max_syn_freq < site.high_freq_corner:
                print("station %s: freq: %f, syn_dt: %f" %
                      (stat, site.high_freq_corner, max_syn_freq))

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

            # Create path names and check if their sizes are within bounds
            sim_file_in = os.path.join(a_outdir_fas,
                                       "%d.%s.smc8.smooth.fs.col" %
                                       (self.sim_id, stat))
            obs_file_in = os.path.join(a_outdir_fas,
                                       "obs.%s.smc8.smooth.fs.col" %
                                       (stat))
            sim_file_tmp = os.path.join(a_tmpdir, "tmp.fas.sim.txt")
            obs_file_tmp = os.path.join(a_tmpdir, "tmp.fas.obs.txt")
            rewrite_fas_eas_file(sim_file_in, sim_file_tmp)
            rewrite_fas_eas_file(obs_file_in, obs_file_tmp)
            outfile = os.path.join(a_outdir, "%s-%d.fas-resid.txt" %
                                   (self.comp_label, self.sim_id))
            bband_utils.check_path_lengths([obs_file_tmp, sim_file_tmp, outfile],
                                           bband_utils.GP_MAX_FILENAME)

            gen_resid_bin = os.path.join(install.A_GP_BIN_DIR,
                                         "gen_resid_tbl_3comp")
            cmd = ("%s bbp_format=1 " % (gen_resid_bin) +
                   "datafile1=%s simfile1=%s " % (obs_file_tmp,
                                                  sim_file_tmp) +
                   "comp1=fash1 comp2=fash2 comp3=seas " +
                   "eqname=%s mag=%s stat=%s lon=%.4f lat=%.4f " %
                   (self.comp_label, self.src_keys['magnitude'],
                    stat, slon, slat) +
                   "vs30=%d cd=%.2f " % (site.vs30, rrup) +
                   "flo=%f fhi=%f " % (1.0 / min(site.high_freq_corner,
                                                 max_syn_freq),
                                       1.0 / site.low_freq_corner) + 
                   "print_header=%d >> %s 2>> %s" %
                   (print_header_fas, outfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Only need to print header the first time
            if print_header_fas == 1:
                print_header_fas = 0

        # Finished per station processing, now summarize and plot the data
        if os.path.exists(fas_resid_output):
            self.summarize_fas(site_list, a_outdir)

        print("FAS GoF Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 8:
        print("Usage: %s " % (PROG_BASE) +
              "source_file station_list magnitude "
              "comp_label method cut_off sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = FASGof(sys.argv[1], sys.argv[2], sys.argv[3],
                sys.argv[4], sys.argv[5],
                cutoff=int(sys.argv[6]),
                sim_id=int(sys.argv[7]))
    ME.run()
