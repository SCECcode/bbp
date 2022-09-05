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

Generic container for the Broadband Platform modules
"""
from __future__ import division, print_function

# Import Python modules
import os
import shutil

# Import Broadband modules
from genslip import Genslip
from ucrmg import UCrmg
from jbsim import Jbsim
from hfsims import Hfsims
from syn1D import Syn1D
from uc_stitch import UCStitch
from bbtoolbox import BBToolbox
from uc_site import UCSite
from wcc_siteamp import WccSiteamp
from rotd50 import RotD50
from fas import FAS
from obs_seismograms import ObsSeismograms
from copy_seismograms import CopySeismograms
from gen_plots import GenPlots
from gp_gof import GPGof
from fas_gof import FASGof
from sdsu_mogof import SDSUMOGoF
from gmpe_plot import GMPEPlot
from gmpe_comparison import GMPEComparison
from calculate_gmpe import CalculateGMPE
from match import Match
from plot_seis import PlotSeis
from plot_map import Plot_Map
from genhtml import GenHTML
from exsim import ExSim
from csm import CSM
from song_rmg_single_seg import SongRMGSS
from song_rmg_multi_seg import SongRMGMS
from as16 import AS16
from rzz2015 import RZZ2015
from rzz2015_gmpe import RZZ2015GMPE
from rotd100 import RotD100
from anderson_gof import AndersonGOF
from irikura_gen_srf import IrikuraGenSrf
from irikura_hf import IrikuraHF
from seismo_soil import SeismoSoil

class Module(object):
    def __init__(self):
        self.module_name = ""
        self.module_args = []
        self.kw_args = dict()
        self.files_to_stage = []

    def setName(self, name):
        self.module_name = name

    def getName(self):
        return self.module_name

    def addArg(self, arg):
        self.module_args.append(arg)

    def addArgs(self, args):
        for arg in args:
            self.module_args.append(arg)

    def setArgs(self, args):
        self.module_args = []
        self.addArgs(args)

    def addKeywordArg(self, keyword, value):
        self.kw_args[keyword] = value

    def addStageFile(self, file_to_stage):
        self.files_to_stage.append(file_to_stage)

    def addStageFiles(self, files):
        for file_to_stage in files:
            self.files_to_stage.append(file_to_stage)

    def resetStageFiles(self):
        self.files_to_stage = []

    def getStageFiles(self):
        return self.files_to_stage

    def stage(self, stage_dir):
        for file_to_stage in self.files_to_stage:
            if os.path.dirname(file_to_stage) == stage_dir:
                # File is already there, skip it
                continue
            if os.path.exists(os.path.join(stage_dir,
                                           os.path.basename(file_to_stage))):
                # File is already there, skip it
                continue
            # print("Staging: %s to %s" % (file, stage_dir))
            shutil.copy2(file_to_stage, stage_dir)

    def getArgs(self):
        return self.module_args

    def getKeywordArgs(self):
        return self.kw_args

    def instantiate(self, sim_id):
        print()
        #print(self.module_name)
        #for arg in self.module_args:
        #       print arg
        #for kw_arg in self.kw_args.keys():
        #       print "keyword %s: value %s" % (kw_arg, self.kw_args[kw_arg])
        #       print kw_arg.__class__
        self.kw_args['sim_id'] = sim_id
        #kwargs = {"simID" : sim_id}
        #return globals()[self.module_name](*self.module_args, simID=my_simID)
        return globals()[self.module_name](*self.module_args, **self.kw_args)
