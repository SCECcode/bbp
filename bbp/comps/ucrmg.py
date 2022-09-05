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

Broadband Platform Version of UCSB Rupture Model Generator
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import uc_fault_utils
from install_cfg import InstallCfg
from ucrmg_cfg import UCrmgCfg

class UCrmg(object):
    """
    Implement UC rupture Model Generator as a Broadband component
    """

    def __init__(self, i_r_velmodel, i_r_srcfile,
                 o_r_srffile, i_vmodel_name, sim_id=0):
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = o_r_srffile
        self.vmodel_name = i_vmodel_name

    def run(self):
        """
        Runs the UCSB rupture generator
        """
        print("UCSB Rupture Generator".center(80, '-'))

        #
        # installation directories
        #
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")
        logfile = os.path.join(install.A_OUT_LOG_DIR,
                               str(sim_id),
                               "%d.ucrmg.log" % (sim_id))
        #
        # Update the Configuration File
        #
        a_srcfile = os.path.join(a_indir, self.r_srcfile)

        cfg = UCrmgCfg(self.vmodel_name, a_srcfile)

        # Make sure the tmp and output directories exist
        bband_utils.mkdirs([a_tmpdir, a_outdir, a_param_outdir],
                           print_cmd=False)

        # Save cwd; change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir)

        # define location of input velocity model file
        # a_velmodel = os.path.join(a_indir, self.r_velmodel)

        # Create ffsp.inp
        r_ffsp_inp = "ffsp.inp"
        a_ffsp_inp = os.path.join(a_tmpdir, r_ffsp_inp)
        uc_fault_utils.uc_create_ffsp_inp(a_ffsp_inp, a_srcfile,
                                          self.vmodel_name)

        # Copy velocity model to work directory
        shutil.copy2(cfg.A_UC_LF_VELMODEL,
                     os.path.join(a_tmpdir,
                                  os.path.basename(cfg.A_UC_LF_VELMODEL)))

        # Save a copy of the ffsp.inp file
        shutil.copy2(a_ffsp_inp, os.path.join(a_param_outdir,
                                              r_ffsp_inp))

        # Now, run the UCSB rupture generator
        cmd = ("%s >> %s 2>&1" %
               (cfg.A_UC_FFSP_EXE, logfile))
        bband_utils.runprog(cmd)

        # Copy output to indata and outdata directories
        ffsp_output_file = "%s.bst" % (cfg.FFSP_OUTPUT_PREFIX)
        shutil.copy2(ffsp_output_file, os.path.join(a_indir, ffsp_output_file))
        shutil.copy2(ffsp_output_file, os.path.join(a_outdir, ffsp_output_file))

        os.chdir(old_cwd)

        print("UCSB Rupture Generator Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = UCrmg(sys.argv[1], sys.argv[2], sys.argv[3],
               sys.argv[4], int(sys.argv[5]))
