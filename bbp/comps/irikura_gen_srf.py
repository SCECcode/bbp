#!/usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import shutil

# Import Broadband modules
import plot_srf
import bband_utils
from irikura_gen_srf_cfg import IrikuraGenSrfCfg
from install_cfg import InstallCfg

class IrikuraGenSrf(object):
    """
    Implements Arben's gen_srf.csh script in Python
    """
    def __init__(self, i_r_velmodel, i_r_srcfile,
                 o_r_srffile, i_vmodel_name, sim_id=0,
                 **kwargs):
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = o_r_srffile
        self.vmodel_name = i_vmodel_name
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
        This function prepares the parameters for Irikura's gen_srf then calls it
        """

        print("IrikuraGenSrf".center(80, '-'))

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
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir,
                            a_logdir, a_param_outdir])

        # Now, file paths
        self.log = os.path.join(a_logdir, "%d.gen_srf.log" % (sim_id))
        a_srcfiles = [os.path.join(a_indir,
                                   srcfile) for srcfile in self.r_srcfiles]

        # Read src file
        cfg = IrikuraGenSrfCfg(a_srcfiles)

        # Define location of input velocity model and output srf file
        if cfg.num_srcfiles > 1:
            a_srffile = os.path.join(a_tmpdir, self.r_srffile)
            a_final_srffile = os.path.join(a_indir, self.r_srffile)
        else:
            a_srffile = os.path.join(a_indir, self.r_srffile)
        a_velmod = os.path.join(install.A_IN_DATA_DIR, str(sim_id),
                                self.r_velmodel)

        # Run in tmpdir subdir to isolate temp fortran files
        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir)

        # Read parameters from the src(s) file(s)
        # The following parameters should be common to all SRC files
        # So we just read from the first one
        simulation_seed = int(cfg.CFGDICT[0]['seed'])
        dip = cfg.CFGDICT[0]['dip']
        rake = cfg.CFGDICT[0]['rake']
        dlen = cfg.CFGDICT[0]['dlen']
        dwid = cfg.CFGDICT[0]['dwid']
        lon_top_center = cfg.CFGDICT[0]['lon_top_center']
        lat_top_center = cfg.CFGDICT[0]['lat_top_center']
        depth_to_top = cfg.CFGDICT[0]['depth_to_top']
        if cfg.num_srcfiles > 1:
            fault_len = cfg.CFGDICT[0]['max_fault_length']
        else:
            fault_len = cfg.CFGDICT[0]['fault_length']
        fault_width = cfg.CFGDICT[0]['fault_width']
        # Average strike of all SRC files
        strike = 0.0
        for segment in range(cfg.num_srcfiles):
            strike = strike + cfg.CFGDICT[segment]['strike']
        strike = math.ceil(strike / cfg.num_srcfiles)
        # Hypocenter (down_dip is common to all src files)
        hypo_down_dip = cfg.CFGDICT[0]['hypo_down_dip']
        if cfg.num_srcfiles > 1:
            hypo_along_stk = 0.0
            for segment in range(cfg.num_srcfiles):
                current_fault_len = cfg.CFGDICT[segment]['fault_length']
                current_hypo_along_stk = cfg.CFGDICT[segment]['hypo_along_stk']
                if abs(current_hypo_along_stk) <= current_fault_len:
                    # Hypocenter in this segment!
                    hypo_along_stk = hypo_along_stk + (current_fault_len / 2.0) + current_hypo_along_stk
                    break
                else:
                    # Not here yet, just add the total length of this segment
                    hypo_along_stk = hypo_along_stk + current_fault_len
            # Now convert hypo_along_stk so that 0.0 is the middle of the fault
            hypo_along_stk = hypo_along_stk - (fault_len / 2.0)
        else:
            hypo_along_stk = cfg.CFGDICT[0]['hypo_along_stk']
        #
        # Run gen_srf code
        #
        progstring = ("%s >> %s 2>&1 << END\n" %
                      (os.path.join(install.A_IRIKURA_BIN_DIR, cfg.GENSRF),
                       self.log) +
                      "%s\n" % a_srffile +
                      "%f %f %f %f %f\n" %
                      (fault_len, fault_width,
                       strike, dip, rake) +
                      "%f %f %f\n" %
                      (lon_top_center, lat_top_center, depth_to_top) +
                      "%f %f\n" % (dlen, dwid) +
                      "%f %f %f %f\n" %
                      (hypo_along_stk, hypo_down_dip,
                       cfg.DENS, cfg.VS) +
                      "%f\n" % (cfg.DT) +
                      "%d\n" % (simulation_seed) +
                      "%s\n" % (a_velmod) +
                      "%f\n" % (cfg.VEL_RUP_FRAC) +
                      "END")
        bband_utils.runprog(progstring)

        # Save single segment srf file for Irikura-2 codes later
        progstring = "cp %s %s" % (a_srffile,
                                   os.path.join(a_tmpdir,
                                                "single_seg.%s" %
                                                (self.r_srffile)))
        bband_utils.runprog(progstring)

        # Tbe segments.midpoint.txt file is created for multi-segment ruptures
        # and also for the Irikura-2 codes that require it for all ruptures

        # Assign the slip from the planar fault to each segment's SRF file
        a_segs_file = os.path.join(a_tmpdir, "segments.midpoint.txt")
        # Write segments' file
        seg_file = open(a_segs_file, 'w')
        seg_file.write("segm  lon       lat    depth   fleng fwidth shypo zhypo strike dip rake\n")
        seg_file.write("%d\n" % (cfg.num_srcfiles))
        total_length = 0.0
        for segment in range(cfg.num_srcfiles):
            if abs(cfg.CFGDICT[segment]['hypo_along_stk']) <= cfg.CFGDICT[segment]['fault_length']:
                hypo_along_stk = cfg.CFGDICT[segment]['hypo_along_stk']
                hypo_down_dip = cfg.CFGDICT[segment]['hypo_down_dip']
            else:
                hypo_along_stk = 999.0
                hypo_down_dip = 999.0
            seg_file.write("seg%d    %.6f %.6f %.1f %.1f %.1f %.1f %.1f %.1f %d %d %d\n" %
                           (segment + 1,
                            cfg.CFGDICT[segment]['lon_top_center'],
                            cfg.CFGDICT[segment]['lat_top_center'],
                            cfg.CFGDICT[segment]['depth_to_top'],
                            total_length,
                            (total_length + cfg.CFGDICT[segment]['fault_length']),
                            cfg.CFGDICT[segment]['fault_width'],
                            hypo_along_stk, hypo_down_dip,
                            cfg.CFGDICT[segment]['strike'],
                            cfg.CFGDICT[segment]['dip'],
                            cfg.CFGDICT[segment]['rake']))
            total_length = total_length + cfg.CFGDICT[segment]['fault_length']
        seg_file.close()

        if cfg.num_srcfiles > 1:
            #
            # Run gen_srf_segment code
            #
            for segment in range(cfg.num_srcfiles):
                progstring = ("%s >> %s 2>&1 << END\n" %
                              (os.path.join(install.A_IRIKURA_BIN_DIR,
                                            cfg.GENSRFSEGMENT), self.log) +
                              ".\n" +
                              "%s\n" % (self.r_srffile) +
                              "./segments.midpoint.txt\n" +
                              "%d\n" % (segment + 1) +
                              "%f %f\n" % (dlen, dwid) +
                              "END")

                # Run code
                bband_utils.runprog(progstring)

            #
            # Now add the segments together
            #
            progstring = ("%s >> %s 2>&1 << END\n" %
                          (os.path.join(install.A_IRIKURA_BIN_DIR,
                                        cfg.SUMSEG), self.log) +
                          ".\n" +
                          "%s\n" % (self.r_srffile) +
                          "./segments.midpoint.txt\n" +
                          "%d\n" % (cfg.num_srcfiles) +
                          "%f %f\n" % (dlen, dwid) +
                          "END")

            # Run code
            bband_utils.runprog(progstring)

            # Copy file to final location
            progstring = "cp %s %s" % (os.path.join(a_tmpdir,
                                                    "all_seg.%s" %
                                                    (self.r_srffile)),
                                       a_final_srffile)
            bband_utils.runprog(progstring)

            # Use copied file from now on
            a_srffile = a_final_srffile

        # Restore working directory
        os.chdir(old_cwd)

        #
        # Move results to output file
        #
        progstring = "cp %s %s" % (a_srffile,
                                   os.path.join(a_tmpdir, self.r_srffile))
        bband_utils.runprog(progstring)
        progstring = "cp %s %s" % (a_srffile,
                                   os.path.join(a_outdir, self.r_srffile))
        bband_utils.runprog(progstring)

        shutil.copy2(os.path.join(a_tmpdir, "stress_drop.out"),
                     os.path.join(a_param_outdir,
                                  "stress_drop.out"))
        shutil.copy2(os.path.join(a_tmpdir, "segments.midpoint.txt"),
                     os.path.join(a_param_outdir,
                                  "segments.midpoint.txt"))


        # Plot SRF
        plot_srf.run(self.r_srffile, sim_id=self.sim_id)

        print("IrikuraGenSrf Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = IrikuraGenSrf(sys.argv[1], sys.argv[2], sys.argv[3],
                       sys.argv[4], sim_id=int(sys.argv[5]))
    ME.run()
