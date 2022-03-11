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

Python module to start the GP rupture generator
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import velocity_models
from install_cfg import InstallCfg
from genslip_cfg import GenslipCfg, calculate_rvfac
import plot_srf

class Genslip(object):
    """
    Implements Genslip as a Broadband component
    Inputs: sim_id, velocity model, and source file
    Outputs: an SRF file
    Assumes that the input directory has been created by the workflow codes,
    and that the files are transfered over into that directory.
    """

    def __init__(self, i_r_velmodel, i_r_srcfile,
                 o_r_srffile, i_vmodel_name, sim_id=0,
                 **kwargs):
        """
        Initialize class variables
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_velmodel = i_r_velmodel
        self.r_srffile = o_r_srffile
        self.vmodel_name = i_vmodel_name
        self.risetime_coef = None
        self.shal_vrup = None
        self.mean_rvfac = None
        self.range_rvfac = None
        self.risetime_fac = None
        self.deep_risetime_fac = None
        self.slip_sigma = None
        self.range_fwidth_frac = None
        self.r_srcfiles = []

        # Get all src files that were passed to us
        if kwargs is not None and len(kwargs) > 0:
            for idx in range(len(kwargs)):
                self.r_srcfiles.append(kwargs['src%d' % (idx)])
        else:
            # Not a multisegment run, just use the single src file
            self.r_srcfiles.append(i_r_srcfile)

    def run(self):
        """
        Runs Genslip
        """
        print("GP Rupture Generator GenSlip".center(80, '-'))

        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir,
                            a_outdir, a_param_outdir], print_cmd=False)

        # Now, file paths
        fault_seg2gsf_bin = os.path.join(install.A_GP_BIN_DIR, "fault_seg2gsf")
        genslip_bin = os.path.join(install.A_GP_BIN_DIR, "genslip-v5.5.2")
        self.log = os.path.join(a_logdir, "%d.genslip.log" % (sim_id))
        a_velfile = os.path.join(a_indir, self.r_velmodel)
        a_srcfiles = [os.path.join(a_indir,
                                   srcfile) for srcfile in self.r_srcfiles]

        # Read source file(s)
        cfg = GenslipCfg(a_srcfiles)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        # Check for velocity model-specific parameters
        vmodel_params = vel_obj.get_codebase_params('gp')
        # Look for RANGE_FWIDTH_FRAC
        if 'RANGE_FWIDTH_FRAC' in vmodel_params:
            self.range_fwidth_frac = float(vmodel_params['RANGE_FWIDTH_FRAC'])
        else:
            self.range_fwidth_frac = cfg.RANGE_FWIDTH_FRAC
        # Look for RISETIME_COEF
        if 'RISETIME_COEF' in vmodel_params:
            self.risetime_coef = float(vmodel_params['RISETIME_COEF'])
        else:
            self.risetime_coef = cfg.RISETIME_COEF
        # Look for SHAL_VRUP
        if 'SHAL_VRUP' in vmodel_params:
            self.shal_vrup = float(vmodel_params['SHAL_VRUP'])
        else:
            self.shal_vrup = cfg.SHAL_VRUP
        # Look for MEAN_RVFAC
        if 'MEAN_RVFAC' in vmodel_params:
            self.mean_rvfac = float(vmodel_params['MEAN_RVFAC'])
        else:
            self.mean_rvfac = cfg.MEAN_RVFAC
        # Look for RANGE_RVFAC
        if 'RANGE_RVFAC' in vmodel_params:
            self.range_rvfac = float(vmodel_params['RANGE_RVFAC'])
        else:
            self.range_rvfac = cfg.RANGE_RVFAC
        # Look for RISETIME_FAC
        if 'RISETIME_FAC' in vmodel_params:
            self.risetime_fac = float(vmodel_params['RISETIME_FAC'])
        else:
            self.risetime_fac = cfg.RISETIME_FAC
        # Look for DEEP_RISETIME_FAC
        if 'DEEP_RISETIME_FAC' in vmodel_params:
            self.deep_risetime_fac = float(vmodel_params['DEEP_RISETIME_FAC'])
        else:
            self.deep_risetime_fac = cfg.DEEP_RISETIME_FAC
        # Look for SLIP SIGMA
        if 'SLIP_SIGMA' in vmodel_params:
            self.slip_sigma = float(vmodel_params['SLIP_SIGMA'])
        else:
            self.slip_sigma = cfg.SLIP_SIGMA

        # Look for DT
        if 'GF_DT' in vmodel_params:
            gf_dt = float(vmodel_params['GF_DT'])
        else:
            raise bband_utils.ParameterError("Cannot find GF_DT parameter in"
                                             "velocity model %s!" %
                                             (self.vmodel_name))

        # If we have multiple SRC files, make sure they are in the right order
        if cfg.num_srcfiles > 1:
            new_dict = []
            for cur_segno in range(1, cfg.num_srcfiles + 1):
                for segment in cfg.CFGDICT:
                    if int(segment["draping_segno"]) == cur_segno:
                        new_dict.append(segment)
                        break
            if len(new_dict) != cfg.num_srcfiles:
                raise bband_utils.ParameterError("Cannot find all segments in"
                                                 "multi-segment event!")
            cfg.CFGDICT = new_dict

        # Calculate nstk,ndip
        mean_fwidth = cfg.CFGDICT[0]["fault_width"]
        range_fwidth = mean_fwidth * self.range_fwidth_frac
        if "common_seed" in cfg.CFGDICT[0]:
            fwidth = calculate_rvfac(mean_fwidth, range_fwidth,
                                     cfg.CFGDICT[0]["common_seed"],
                                     count=1)
        else:
            fwidth = calculate_rvfac(mean_fwidth, range_fwidth,
                                     cfg.CFGDICT[0]["seed"],
                                     count=1)
        flen_total = 0.0
        for segment in cfg.CFGDICT:
            flen_total = flen_total + segment["fault_length"]
        nstk_tot = round(flen_total / cfg.CFGDICT[0]["dlen"])
        ndip = int(round(fwidth / cfg.CFGDICT[0]["dwid"]))
        fwidth = ndip * cfg.CFGDICT[0]["dwid"]

        # Check if we are using the draping method
        if ((not "draping_group" in cfg.CFGDICT[0]) or
            ("draping_group" in cfg.CFGDICT[0] and
             cfg.CFGDICT[0]["draping_group"] < 1)):
            # For single-segmenmt or multi-segment SRC files
            # not using the draping method

            # Set Hypocenter
            master_shypo = cfg.CFGDICT[0]["hypo_along_stk"]

            # Calculate rvfac
            if "common_seed" in cfg.CFGDICT[0]:
                rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                        cfg.CFGDICT[0]["common_seed"])
            else:
                rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                        cfg.CFGDICT[0]["seed"])

            if "rupture_delay" in cfg.CFGDICT[0]:
                rupture_delay = cfg.CFGDICT[0]["rupture_delay"]
            else:
                rupture_delay = 0.0

            if "moment_fraction" in cfg.CFGDICT[0]:
                moment_fraction = cfg.CFGDICT[0]["moment_fraction"]
            else:
                moment_fraction = -1.0

            if "max_fault_length" in cfg.CFGDICT[0]:
                flen_max = cfg.CFGDICT[0]["max_fault_length"]
            else:
                flen_max = -1.0
        else:
            rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                    cfg.CFGDICT[0]["seed"])
            # Disable the following parameters when using draping
            rupture_delay = 0.0
            moment_fraction = -1.0
            flen_max = -1.0

            # Now set hypocenter along strike
            master_shypo = 0.0
            for segment in cfg.CFGDICT:
                if segment["true_hypo"] != 1:
                    master_shypo = master_shypo + segment["fault_length"]
                else:
                    master_shypo = (master_shypo +
                                    0.5 * (segment["fault_length"] - flen_total) +
                                    segment["hypo_along_stk"])
                    break

        r_gsftmp = "m%.2f-%.2fx%.2f.gsf" % (cfg.CFGDICT[0]["magnitude"],
                                            cfg.CFGDICT[0]["dlen"],
                                            cfg.CFGDICT[0]["dwid"])
        a_fault_seg_in = os.path.join(a_tmpdir, "fault_seg.in")
        a_gsftmp = os.path.join(a_tmpdir, r_gsftmp)

        r_outroot = "m%.2f-%.2fx%.2f_s%d-v5.4.1" % (cfg.CFGDICT[0]["magnitude"],
                                                    cfg.CFGDICT[0]["dlen"],
                                                    cfg.CFGDICT[0]["dwid"],
                                                    cfg.CFGDICT[0]["seed"])
        a_srffile = os.path.join(a_indir, "%s.srf" % (r_outroot))

        # Write fault_seg.in file
        fault_seg = open(a_fault_seg_in, 'w')
        fault_seg.write("%d\n" % (cfg.num_srcfiles))
        for segment in cfg.CFGDICT:
            nstk_seg = int(round(segment["fault_length"] / segment["dlen"]))
            fault_seg.write("%f %f %f %f %f %f %f %.1f %d %d\n" %
                            (segment["lon_top_center"],
                             segment["lat_top_center"],
                             segment["depth_to_top"],
                             segment["strike"],
                             segment["dip"],
                             segment["rake"],
                             segment["fault_length"],
                             fwidth,
                             nstk_seg,
                             ndip))
        fault_seg.close()

        # Save fault_seg_in file
        shutil.copy2(a_fault_seg_in,
                     os.path.join(a_param_outdir,
                                  os.path.basename(a_fault_seg_in)))

        # Output parameters used for this run
        print("Sim parameters: id: %d - fwidth: %.1f - rvfrac: %f" %
              (sim_id, fwidth, rvfac))

        progstring = ("%s read_slip_vals=0 < %s > %s 2>> %s\n" %
                      (fault_seg2gsf_bin, a_fault_seg_in, a_gsftmp, self.log))
        bband_utils.runprog(progstring)

        progstring = ("%s read_erf=0 write_srf=1 " %
                      (genslip_bin) +
                      "read_gsf=1 write_gsf=0 infile=%s " % (a_gsftmp) +
                      "mag=%.2f nstk=%d ndip=%d " %
                      (cfg.CFGDICT[0]["magnitude"], nstk_tot, ndip) +
                      "ns=1 nh=1 " +
                      "kmodel=2 seed=%d slip_sigma=%f " %
                      (cfg.CFGDICT[0]["seed"], self.slip_sigma) +
                      "circular_average=0 modified_corners=0 " +
                      "velfile=%s shypo=%f dhypo=%f rvfrac=%f " %
                      (a_velfile, master_shypo,
                       cfg.CFGDICT[0]["hypo_down_dip"], rvfac) +
                      "shal_vrup_dep=%f shal_vrup_deprange=%f shal_vrup=%f " %
                      (cfg.RTDEP, cfg.RTDEP_RANGE, self.shal_vrup) +
                      "side_taper=0.02 bot_taper=0.0 top_taper=0.0 " +
                      "dt=%f risetime_coef=%f plane_header=1 " %
                      (gf_dt, self.risetime_coef) +
                      "risetimefac=%f risetimedep=%f risetimedep_range=%f " %
                      (self.risetime_fac, cfg.RTDEP, cfg.RTDEP_RANGE) +
                      "rt_scalefac=%f slip_water_level=%f " %
                      (cfg.RT_SCALEFAC, cfg.SLIP_WATER_LEVEL) +
                      "deep_risetimedep=%f deep_risetimedep_range=%f " %
                      (cfg.DEEP_RISETIMEDEP, cfg.DEEP_RISETIMEDEP_RANGE) +
                      "deep_risetimefac=%f " %
                      (self.deep_risetime_fac) +
                      "flen_max=%f rupture_delay=%f moment_fraction=%f " %
                      (flen_max, rupture_delay, moment_fraction) +
                      "srf_version=2.0 rake_sigma=15.0 fdrup_time=1 " +
                      "deep_vrup=0.6 use_gaus=1 alpha_rough=0.01 " +
                      "lambda_min=0.08 tsfac_coef=1.1 tsfac1_sigma=1.0 " +
                      "tsfac1_scor=0.8 rtime1_sigma=%f rtime1_scor=0.5 " %
                      (self.slip_sigma) +
                      "tsfac_bzero=-0.1 tsfac_slope=-0.5 " +
                      "> %s 2>> %s" % (a_srffile, self.log))
        bband_utils.runprog(progstring)

        #
        # mv result to outputfile
        #
        progstring = "cp %s %s" % (a_srffile, os.path.join(a_tmpdir, self.r_srffile))
        bband_utils.runprog(progstring)
        progstring = "cp %s %s" % (a_srffile, os.path.join(a_indir, self.r_srffile))
        bband_utils.runprog(progstring)
        progstring = "cp %s %s" % (a_srffile, os.path.join(a_outdir, self.r_srffile))
        bband_utils.runprog(progstring)

        # Plot SRF
        plot_srf.run(self.r_srffile, sim_id=self.sim_id)

        print("GP GenSlip Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = Genslip(sys.argv[1], sys.argv[2],
                 sys.argv[3], sys.argv[4],
                 sim_id=int(sys.argv[5]))
    sys.exit(0)
