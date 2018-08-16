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

Broadband Platform Version of Martin Mai BBcoda2.csh
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
# import math # (used to calculate moment)

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
                 o_r_srffile, i_vmodel_name, sim_id=0):
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

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir], print_cmd=False)

        # Now, file paths
        self.log = os.path.join(a_logdir, "%d.genslip.log" % (sim_id))
        a_srcfile = os.path.join(a_indir, self.r_srcfile)
        a_velfile = os.path.join(a_indir, self.r_velmodel)

        # Read src file
        cfg = GenslipCfg(a_srcfile)

        # Define location of input velocity model file
        a_velmodel = os.path.join(a_tmpdir, self.r_velmodel)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        # Check for velocity model-specific parameters
        vmodel_params = vel_obj.get_codebase_params('gp')
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

        # Calculate nstk,ndip
        nstk = round(cfg.CFGDICT["fault_length"] / cfg.CFGDICT["dlen"])
        ndip = round(cfg.CFGDICT["fault_width"] / cfg.CFGDICT["dwid"])

        # Calculate rvfac
        if "common_seed" in cfg.CFGDICT:
            rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                    cfg.CFGDICT["common_seed"])
        else:
            rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                    cfg.CFGDICT["seed"])

        # moment = math.pow(10, 1.5 * (cfg.MAG + 10.7))

        # For multi-segment SRC files
        if "rupture_delay" in cfg.CFGDICT:
            rupture_delay = cfg.CFGDICT["rupture_delay"]
        else:
            rupture_delay = 0.0

        if "moment_fraction" in cfg.CFGDICT:
            moment_fraction = cfg.CFGDICT["moment_fraction"]
        else:
            moment_fraction = -1.0

        if "max_fault_length" in cfg.CFGDICT:
            flen_max = cfg.CFGDICT["max_fault_length"]
        else:
            flen_max = -1.0

        r_gsftmp = "m%.2f-%.2fx%.2f.gsf" % (cfg.CFGDICT["magnitude"],
                                            cfg.CFGDICT["dlen"],
                                            cfg.CFGDICT["dwid"])
        a_gsftmp = os.path.join(a_tmpdir, r_gsftmp)

        r_outroot = "m%.2f-%.2fx%.2f_s%d-v5.2.2" % (cfg.CFGDICT["magnitude"],
                                                    cfg.CFGDICT["dlen"],
                                                    cfg.CFGDICT["dwid"],
                                                    cfg.CFGDICT["seed"])
        a_srffile = os.path.join(a_indir, "%s.srf" % (r_outroot))

        progstring = ("%s/fault_seg2gsf read_slip_vals=0 << EOF > %s 2>> %s\n" %
                      (install.A_GP_BIN_DIR, a_gsftmp, self.log) +
                      "1\n" +
                      "%f %f %f %f %f %f %f %f %d %d\n" %
                      (cfg.CFGDICT["lon_top_center"],
                       cfg.CFGDICT["lat_top_center"],
                       cfg.CFGDICT["depth_to_top"],
                       cfg.CFGDICT["strike"], cfg.CFGDICT["dip"],
                       cfg.CFGDICT["rake"], cfg.CFGDICT["fault_length"],
                       cfg.CFGDICT["fault_width"], nstk, ndip) + "EOF")
        bband_utils.runprog(progstring)

        progstring = ("%s/genslip-v5.2.2 read_erf=0 write_srf=1 " %
                      (install.A_GP_BIN_DIR) +
                      "read_gsf=1 write_gsf=0 infile=%s " % (a_gsftmp) +
                      "mag=%.2f nstk=%d ndip=%d " %
                      (cfg.CFGDICT["magnitude"], nstk, ndip) +
                      "ns=1 nh=1 " +
                      "kmodel=2 seed=%d slip_sigma=%f " %
                      (cfg.CFGDICT["seed"], self.slip_sigma) +
                      "circular_average=0 modified_corners=0 " +
                      "velfile=%s shypo=%f dhypo=%f rvfrac=%f " %
                      (a_velfile, cfg.CFGDICT["hypo_along_stk"],
                       cfg.CFGDICT["hypo_down_dip"], rvfac) +
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
                      "tsfac1_scor=0.8 rtime1_sigma=0.85 rtime1_scor=0.8 " +
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
