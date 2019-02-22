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

Broadband Platform Version of Rob hfsim-stats
Outputs velocity (cm/s)
"""
from __future__ import division, print_function

# Import Python modules
import ast
import os
import sys

# Import Broadband modules
import bband_utils
import validation_cfg
import velocity_models
from hfsims_cfg import HfsimsCfg, calculate_rvfac
from install_cfg import InstallCfg
from station_list import StationList

class Hfsims(object):
    """
    Implement Robert Graves hfsim.csh as a python component
    """
    def __init__(self, i_r_velmodel, i_r_srcfile, i_r_srffile, i_r_stations,
                 i_vmodel_name, val_name=None, sim_id=0):
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_stations = i_r_stations
        self.vmodel_name = i_vmodel_name
        self.val_name = val_name
        self.val_obj = None
        self.qfexp = None
        self.kappa = None
        self.c0 = None
        self.c1 = None
        self.sdrop = None
        self.default_dx = None
        self.default_dy = None
        self.default_fcfac = None
        self.rayset = None
        self.tlen = None
        self.dt = None
        self.mean_rvfac = None
        self.range_rvfac = None
        self.shal_rvfac = None
        self.vsmoho = None
        self.path_dur_model = None
        self.deep_rvfac = None
        self.rvsig = None

    def run(self):
        """
        This function prepares the parameters for HFSim and then calls it
        """
        print("GP HfSims".center(80, '-'))

        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Find validation object if this is a validation run
        if self.val_name is not None:
            self.val_obj = validation_cfg.VE_EVENTS.get_event_by_name(self.val_name)

        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.hfsims_%s.log" % (sim_id, sta_base))

        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))

        a_srffile = os.path.join(a_indir, self.r_srffile)
        # Make sure we work when starting from an SRF file
        if self.r_srcfile:
            a_srcfile = os.path.join(a_indir, self.r_srcfile)
        else:
            a_srcfile = ""
        # Set up basic parameters, read SRC file if provided
        config = HfsimsCfg(a_srcfile=a_srcfile)

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))

        # Check for velocity model-specific parameters
        vmodel_params = vel_obj.get_codebase_params('gp')
        # Look for KAPPA
        if 'KAPPA' in vmodel_params:
            self.kappa = float(vmodel_params['KAPPA'])
        else:
            self.kappa = config.KAPPA
        # Look for QFEXP
        if 'QFEXP' in vmodel_params:
            self.qfexp = float(vmodel_params['QFEXP'])
        else:
            self.qfexp = config.DEFAULT_QFEXP
        # Look for SDROP
        if 'SDROP' in vmodel_params:
            self.sdrop = int(vmodel_params['SDROP'])
        else:
            self.sdrop = config.DEFAULT_SDROP
        # Look for C0 and C1
        if 'C0' in vmodel_params:
            self.c0 = int(vmodel_params['C0'])
        else:
            self.c0 = config.DEFAULT_C0
        if 'C1' in vmodel_params:
            self.c1 = int(vmodel_params['C1'])
        else:
            self.c1 = config.DEFAULT_C1
        # Look for DEFAULT_FCFAC
        if 'DEFAULT_FCFAC' in vmodel_params:
            self.default_fcfac = float(vmodel_params['DEFAULT_FCFAC'])
        else:
            self.default_fcfac = config.DEFAULT_FCFAC
        # Look for rayset
        if 'RAYSET' in vmodel_params:
            self.rayset = ast.literal_eval(vmodel_params['RAYSET'])
        else:
            self.rayset = config.RAYSET
        # Look for a high frequency DT
        if 'HF_DT' in vmodel_params:
            self.dt = float(vmodel_params['HF_DT'])
        else:
            self.dt = config.DT
        # Look for MEAN_RVFAC
        if 'MEAN_RVFAC' in vmodel_params:
            self.mean_rvfac = float(vmodel_params['MEAN_RVFAC'])
        else:
            self.mean_rvfac = config.MEAN_RVFAC
        # Look for RANGE_RVFAC
        if 'RANGE_RVFAC' in vmodel_params:
            self.range_rvfac = float(vmodel_params['RANGE_RVFAC'])
        else:
            self.range_rvfac = config.RANGE_RVFAC
        # Look for SHAL_RVFAC
        if 'SHAL_RVFAC' in vmodel_params:
            self.shal_rvfac = float(vmodel_params['SHAL_RVFAC'])
        else:
            self.shal_rvfac = config.SHAL_RVFAC
        # Look for VSMOHO
        if 'VSMOHO' in vmodel_params:
            self.vsmoho = float(vmodel_params['VSMOHO'])
        else:
            self.vsmoho = config.DEFAULT_VSMOHO
        # Look for DEEP_RVFAC
        if 'DEEP_RVFAC' in vmodel_params:
            self.deep_rvfac = float(vmodel_params['DEEP_RVFAC'])
        else:
            self.deep_rvfac = config.DEEP_RVFAC
        # Look for PATH_DUR_MODEL
        if 'PATH_DUR_MODEL' in vmodel_params:
            self.path_dur_model = int(vmodel_params['PATH_DUR_MODEL'])
        else:
            self.path_dur_model = config.PATH_DUR_MODEL
        # Look for RVSIG
        if 'RVSIG' in vmodel_params:
            self.rvsig = float(vmodel_params['RVSIG'])
        else:
            self.rvsig = config.RVSIG
        # Look for DX
        if 'DEFAULT_DX' in vmodel_params:
            self.default_dx = float(vmodel_params['DEFAULT_DX'])
        else:
            self.default_dx = config.DEFAULT_DX
        # Look for DY
        if 'DEFAULT_DY' in vmodel_params:
            self.default_dy = float(vmodel_params['DEFAULT_DY'])
        else:
            self.default_dy = config.DEFAULT_DY
        # Look for ISPAR_ADJUST
        if 'ISPAR_ADJUST' in vmodel_params:
            ispar_adjust = int(vmodel_params['ISPAR_ADJUST'])
        else:
            ispar_adjust = config.ISPAR_ADJUST

        # Calculate rvfac
        if "common_seed" in config.CFGDICT:
            rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                    config.CFGDICT["common_seed"])
        else:
            rvfac = calculate_rvfac(self.mean_rvfac, self.range_rvfac,
                                    config.CFGDICT["seed"])

        # Look for tlen
        if "TLEN" in vmodel_params:
            self.tlen = float(vmodel_params['TLEN'])
        else:
            self.tlen = config.TLEN

        # Start with some default values
        moment = -1
        extra_fcfac = config.DEFAULT_EXTRA_FCFAC

        if self.val_obj is not None:
            extra_fcfac = float(self.val_obj.get_input("GP", "EXTRA_FCFAC"))
            try:
                tlen = float(self.val_obj.get_input("GP", "TLEN"))
                self.tlen = tlen
            except (ValueError, KeyError, TypeError):
                # No problem, just use the default TLEN for this simulation
                pass

        fcfac = round((1 + self.default_fcfac) * (1 + extra_fcfac) - 1, 4)

        a_slipfile = os.path.join(a_tmpdir, "%s.%s.%fx%f" % (self.r_srffile,
                                                             sta_base,
                                                             self.default_dx,
                                                             self.default_dy))

        progstring = ("%s " %
                      (os.path.join(install.A_GP_BIN_DIR, "srf2stoch")) +
                      "infile=%s outfile=%s " %
                      (a_srffile, a_slipfile) +
                      "target_dx=%f target_dy=%f " %
                      (self.default_dx, self.default_dy) +
                      ">> %s 2>&1" % (self.log))
        bband_utils.runprog(progstring)

        a_outp = os.path.join(a_tmpdir, "tmp_hfsim_out")
        a_velmod = os.path.join(install.A_IN_DATA_DIR, str(sim_id),
                                self.r_velmodel)
        a_statfile = os.path.join(install.A_IN_DATA_DIR, str(sim_id),
                                  self.r_stations)

        # Create local velocity model
        vel_in_fp = open(a_velmod, 'r')
        a_velmod = "%s_%s.local" % (a_velmod, sta_base)
        vel_out_fp = open(a_velmod, 'w')
        vel_in_data = vel_in_fp.readlines()
        vel_in_fp.close()
        i = 0
        for line in vel_in_data:
            i += 1
            if line.startswith('#') or line.startswith('%'):
                continue
            pieces = line.split()
            if len(pieces) >= 4:
                th = float(pieces[0])
                vp = float(pieces[1])
                vs = float(pieces[2])
                dn = float(pieces[3])
                qs = self.c0 + self.c1 * vs
                if i == len(vel_in_data):
                    th = 0.0
                vel_out_fp.write("%8.4f %8.4f %8.4f %8.4f %8.2f %8.2f\n" %
                                 (th, vp, vs, dn, qs, qs))
            else:
                vel_out_fp.write(line)
        vel_out_fp.flush()
        vel_out_fp.close()
        #
        # Scan the station list with this object construction
        # This scanner removes all the comment lines and the
        # list that is returned has one station per line in it.
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()
        nstat = len(site_list)

        # Create rayset param list
        rayset_param = ""
        for item in self.rayset:
            rayset_param = rayset_param + "%d " % (item)
        rayset_param = rayset_param.strip()
        #
        # Run initial hfsim conf
        #
        progstring = ("%s >> %s 2>&1 << END\n" %
                      (os.path.join(install.A_GP_BIN_DIR, config.HFSIM),
                       self.log) +
                      "%d\n" % self.sdrop +
                      "%s\n" % a_statfile +
                      "%s\n" % a_outp +
                      "%s\n" % rayset_param +
                      "%d\n" % config.SITEAMP +
                      "4 0 0.02 19.9\n" +
                      "%d\n" % config.CFGDICT["seed"] +
                      "%d\n" % nstat +
                      "%f %f %f %f %f\n" %
                      (self.tlen, self.dt, config.FMAX,
                       self.kappa, self.qfexp) +
                      "%f %f %f %f %f\n" %
                      (rvfac, self.shal_rvfac, self.deep_rvfac,
                       config.C_ZERO, config.C_ALPHA) +
                      "%s %f\n" % (moment, config.RUPV) +
                      "%s\n" % a_slipfile +
                      "%s\n" % a_velmod +
                      "%f\n" % self.vsmoho +
                      "-99 0.0 0.0 0.0 0.0 1\n" +
                      "-1\n" +
                      "%f 0.0 %f\n" % (config.FA_SIG1, self.rvsig) +
                      "%d\n" % (self.path_dur_model) +
                      "%d -1 -1\n" % (ispar_adjust) +
                      "END")
        bband_utils.runprog(progstring)

        #
        # Start the per station processing
        #
        for site in site_list:
            # Need to integrate each component, since hfsims outputs cm/s/s
            for comp in ['000', '090', 'ver']:
                cmd = ("%s integ=1 " %
                       (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                       "filein=%s_%s.%s fileout=%s/%d.%s-hf.%s >> %s 2>&1" %
                       (a_outp, site.scode, comp, a_tmpdir, sim_id,
                        site.scode, comp, self.log))
                bband_utils.runprog(cmd, print_cmd=False)

            progstring1 = ("%s wcc2bbp=1 " %
                           (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                           'title="HF Sim NGAH, stat=%s" ' % (site.scode) +
                           "nsfile=%s/%d.%s-hf.000 " %
                           (a_tmpdir, sim_id, site.scode) +
                           "ewfile=%s/%d.%s-hf.090 " %
                           (a_tmpdir, sim_id, site.scode) +
                           "udfile=%s/%d.%s-hf.ver " %
                           (a_tmpdir, sim_id, site.scode) +
                           'units="%s" > %s/%d.%s-hf.bbp 2>> %s\n' %
                           (str(config.UNITS), a_tmpdir, sim_id,
                            site.scode, self.log))
            bband_utils.runprog(progstring1, print_cmd=False)

            progstring1 = ("%s wcc2bbp=1 " %
                           (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                           'title="HF Sim NGAH, stat=%s" ' % site.scode +
                           "nsfile=%s_%s.000 " %
                           (a_outp, site.scode) +
                           "ewfile=%s_%s.090 " %
                           (a_outp, site.scode) +
                           "udfile=%s_%s.ver " %
                           (a_outp, site.scode) +
                           '> %s_%s.bbp 2>> %s\n' %
                           (a_outp, site.scode, self.log))
            bband_utils.runprog(progstring1, print_cmd=False)

        print("GP HfSims Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = Hfsims(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
                sys.argv[6], sim_id=int(sys.argv[7]))
    ME.run()
