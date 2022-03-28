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

This modules defines the configuration parameters for the SDSU mogof script
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
from install_cfg import InstallCfg

class SDSUMOGofCfg(object):
    """
    Define the configuration parameters for the SDSU MO-GOF program
    """

    def __init__(self):

        install = InstallCfg.getInstance()
        # Name and Path to executable

        self.COMPS = ["000", "090", "ver"]
        self.INPUT_FORMAT_GOF = "GOF"

        # Path to GOF binaries
        self.GOF_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_MO")
        self.GOF_SpFit_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_SpFit")
        self.GOF_PGX_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_PGX")
        self.GOF_CCFit_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_CCFit")
        self.GOF_DCumEn_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_DCumEn")
        self.GOF_FSComp_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_FSComp")
        self.GOF_InElFit_BIN = os.path.join(install.A_SDSU_BIN_DIR,
                                            "GOF_InElFit")
        self.GOF_SAFit16_BIN = os.path.join(install.A_SDSU_BIN_DIR,
                                            "GOF_SAFit16")
        self.GOF_SpecDurFit_BIN = os.path.join(install.A_SDSU_BIN_DIR,
                                               "GOF_SpecDurFit")
        self.GOF_NGA_BIN = os.path.join(install.A_SDSU_BIN_DIR, "GOF_MO_NGA")

        # Working file names
        self.PARAM_DAT_FILE = "PARAM.dat"
        self.INPUT_SEISMO_1 = "CONCAT_1"
        self.INPUT_SEISMO_2 = "CONCAT_2"

        # Summary file
        self.SUMMARY_FILE = "gof_summary.txt"

        self.cfggof = {}
        self.cfggof["input_format"] = "GOF"
        self.cfggof["input_set_1"] = "./sample1_1.input"
        self.cfggof["input_set_2"] = "./sample1_2.input"
        self.cfggof["output_dir"] = "./output"
        self.cfggof["work_dir"] = "./work"
        self.cfggof["use_nga"] = False

        self.cfggof["num_headers"] = 0 #number of headers preceeding each seismogram
        self.cfggof["num_station"] = 1 #number of stations
        self.cfggof["timesteps_set_1"] = 0 #nt number of timesteps in set 1
        self.cfggof["timesteps_set_2"] = 0 #nt number of timesteps in set 2
        self.cfggof["input_param"] = "A" #Input format(A,V,D), (meter, sec)
        self.cfggof["seismo_length"] = 0 #length of seismograms (seconds)
        self.cfggof["low_cut"] = 0.1 #Low cut (period)
        self.cfggof["high_cut"] = 10.0 #High cut (period)
        self.cfggof["weights"] = dict()
        self.cfggof["weights"]["pga"] = 1.0                          #  Weighting on PGA
        self.cfggof["weights"]["pgv"] = 1.0                          #  Weighting on PGV
        self.cfggof["weights"]["pgd"] = 1.0                          #  Weighting on PGD
        self.cfggof["weights"]["psa"] = 1.0                          #  Weighting on PSA
        self.cfggof["weights"]["spectral_Fit"] = 1.0                 #  Weighting on Spectral Fit
        self.cfggof["weights"]["cumulative_energy_fit"] = 1.0        #  Weighting on Cumulative Energy Fit
        self.cfggof["weights"]["inelastic_elastic_fit"] = 1.0        #  Weighting on Inelastic/Elastic Fit (16)
        self.cfggof["weights"]["sepctral_acc"] = 1.0                 #  Weighting on Spec Acc (16)
        self.cfggof["weights"]["spec_duration"] = 1.0                #  Weighting on Spec Dur (16)
        self.cfggof["weights"]["data_energy_release_duration"] = 1.0 #  Weighting on Data Energy Release Duration (5%-75%)
        self.cfggof["weights"]["cross_correlation"] = 1.0            #  Weighting on Cross-Correlation
        self.cfggof["weights"]["fourier_spectrum"] = 1.0             #  Weighting on Fourier Spectrum
        self.cfggof["file"] = dict()
        self.cfggof["file"]["pga"] = "GOF_PGA.list"
        self.cfggof["file"]["pgv"] = "GOF_PGV.list"
        self.cfggof["file"]["pgd"] = "GOF_PGD.list"
        self.cfggof["file"]["psa"] = "GOF_PSA.list"
        self.cfggof["file"]["spectral_Fit"] = "GOF_SPECFIT.list"
        self.cfggof["file"]["cumulative_energy_fit"] = "GOF_ENERGYFIT.list"
        self.cfggof["file"]["inelastic_elastic_fit"] = "GOF_InElEl.list"
        self.cfggof["file"]["sepctral_acc"] = "GOF_SAFIT.list"
        self.cfggof["file"]["spec_duration"] = "GOF_SPECDUR.list"
        self.cfggof["file"]["data_energy_release_duration"] = "GOF_DUR.list"
        self.cfggof["file"]["cross_correlation"] = "GOF_CROSSCOR.list"
        self.cfggof["file"]["fourier_spectrum"] = "GOF_FS.list"

        self.cfgnga = {}
        self.cfgnga["source_mag"] = 0.0
        self.cfgnga["dip"] = 0.0
        self.cfgnga["rake"] = 0.0
        self.cfgnga["depth_coseismic"] = 0.0
        self.cfgnga["site_file"] = ""

if __name__ == "__main__":
    print("Test Config Class: %s" % (sys.argv[0]))
