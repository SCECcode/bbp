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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import fault_utils
import stas2files
from install_cfg import InstallCfg
from station_list import StationList

class UCStitch(object):
    '''
    To be run after using syn1D for HF and someone else for LF
    Uses Jan's wavelet algorithm to combine HF and LF data
    '''

    def __init__(self, i_r_velocity, i_r_srcfile, i_r_srffile,
                 i_r_stations, acc=False, sim_id=0):
        self.sim_id = sim_id
        self.r_velocity = i_r_velocity
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_stations = i_r_stations
        self.acc = acc

    def make_bbp(self, a_tmpdir, site):
        print("Making bbp for %s." % (site))
        bbp_fp = open(os.path.join(a_tmpdir,
                                   "%d.%s-stitch.bbp" % (self.sim_id, site)),
                      'w')

        gm_filenames = []
        gm_filenames.append("%s/%s.000.gmBB.001" % (a_tmpdir, site))
        gm_filenames.append("%s/%s.090.gmBB.001" % (a_tmpdir, site))
        gm_filenames.append("%s/%s.ver.gmBB.001" % (a_tmpdir, site))
        gm_data = []
        for gm_filename in gm_filenames:
            print(gm_filename)
            tmp_fp = open(gm_filename, 'r')
            tmp_data = []
            header = tmp_fp.readline()
            for line in tmp_fp.readlines():
                pieces = line.split()
                for piece in pieces:
                    tmp_data.append(piece.strip())
            gm_data.append(tmp_data)
            tmp_fp.close()

        nt = int(header.split()[0])
        dt = float(header.split()[1])
        ts = 0.0
        print(nt, dt)

        bbp_fp.write("# nt = %d   dt = %f\n" % (nt, dt))
        bbp_fp.write("# t\tn/s\te/w\tu/d\n")

        for i in range(0, len(gm_data[0])):
            bbp_fp.write("%e\t%f\t%f\t%f\n" %
                         (ts, float(gm_data[0][i]),
                          float(gm_data[1][i]), float(gm_data[2][i])))
            ts += dt
        bbp_fp.flush()
        bbp_fp.close()

    def split_bbp(self, a_tmpdir, site, bbp_file):
        #splits BBP file into nl1D file
        print(bbp_file)
        bbp_fp = open(bbp_file, 'r')
        bbp_data = bbp_fp.readlines()
        bbp_fp.close()

        nl_filenames = []
        nl_filenames.append("%s/%s.000.nl1D.001" % (a_tmpdir, site))
        nl_filenames.append("%s/%s.090.nl1D.001" % (a_tmpdir, site))
        nl_filenames.append("%s/%s.ver.nl1D.001" % (a_tmpdir, site))
        nl_fp = []
        for nl_filename in nl_filenames:
            nl_fp.append(open(nl_filename, 'w'))
        line = bbp_data[0]
        while line.startswith("%") or line.startswith("#"):
            del bbp_data[0]
            line = bbp_data[0]
        for fp in nl_fp:
            #nt dt
            fp.write("\t%d\t%f\n" % (len(bbp_data),
                                     float(bbp_data[1].split()[0])))
        count = 0
        for line in bbp_data:
            pieces = line.split()
            for i in range(0, 3):
                nl_fp[i].write("%f " % float(pieces[i+1]))
            count += 1
            if (count % 5) == 0:
                for i in range(0, 3):
                    nl_fp[i].write("\n")
        for fp in nl_fp:
            fp.flush()
            fp.close()

    def run(self):
        """
        Run the UCSB Stitch code
        """
        print("UCSB Stitch".center(80, '-'))

        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.uc_stitch_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR,
                                    str(self.sim_id),
                                    "uc_stitch_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        #
        # Make sure the outpute and tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir],
                           print_cmd=False)

        a_velocity = os.path.join(a_indir, self.r_velocity)
        print("UC_stich - SRF File: %s" % (self.r_srffile))
        if not os.path.isfile(self.r_srffile):
            a_srffile = os.path.join(a_indir, self.r_srffile)
            if not os.path.isfile(a_srffile):
                print("Error (uc_stich): Unable to locate SRF file:",
                      a_srffile,
                      "\nExiting Broadband...")
                sys.exit(1)
        else:
            a_srffile = self.r_srffile

        # Copy srf and velocity files to tmpdir
        shutil.copy2(a_srffile, os.path.join(a_tmpdir_mod, self.r_srffile))
        shutil.copy2(a_velocity, os.path.join(a_tmpdir_mod, self.r_velocity))

        # Store cwd and change over to tmpdir so the executable can
        # find the files
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        a_stations = os.path.join(a_indir, self.r_stations)
        print(a_stations)
        slo = StationList(a_stations)
        site_list = slo.get_station_list()

        # Need stations in UCSB LL format
        a_uc_stations = os.path.join(a_tmpdir_mod, "stations.ll")
        a_uc_vs30 = os.path.join(a_tmpdir_mod, "stations.vs30")
        stas2files.gp2uc_stalist(slo, a_uc_stations, a_uc_vs30)

        # Convert stations to UCSB XY format, if it doesn't exist
        if not os.path.exists("stations.xy"):
            r_faultfile = "faultGlobal.in"
            a_faultfile = os.path.join(a_tmpdir_mod, r_faultfile)
            if not os.path.isfile(a_faultfile):
                if os.path.isfile(os.path.join(a_indir, r_faultfile)):
                    shutil.copy2(os.path.join(a_indir, r_faultfile),
                                 a_faultfile)
                else:
                    print("Extracting faultGlobal.in from "
                          "SRF and velocity model.")
                    fp = open("faultGlobalTmp", 'w')
                    fp.write("%s\n" % self.r_velocity)
                    fp.write("%s\n" % self.r_srffile)
                    fp.close()
                    cmd = ("%s/getfaultGlobal < faultGlobalTmp >> %s 2>&1" %
                           (install.A_UCSB_BIN_DIR, self.log))
                    bband_utils.runprog(cmd)

            cmd = "%s/statLL2XY >> %s 2>&1 " % (install.A_UCSB_BIN_DIR,
                                                self.log)
            bband_utils.runprog(cmd)

        #write stitch BBP
        sBBP_in = open("stitchBBP.inp", "w")
        #stations list
        sBBP_in.write("stations.xy\n")
        #velocity file list
        sBBP_in.write("VMname.list\n")
        #dir with 1D synthetics (<stat>.3comp)
        sBBP_in.write("%s/\n" % a_tmpdir)
        #dir with non-linear 3D synthetics (<stat>.<comp>.nl1D.001)
        sBBP_in.write("%s/\n" % a_tmpdir)
        #joint frequency, fmax
        sBBP_in.write("1.0, 15.0\n")
        # depth of hypocenter
        if self.r_srcfile is not None and self.r_srcfile != "":
            a_srcfile = os.path.join(a_indir, self.r_srcfile)
            hypo_dep = fault_utils.calculate_hypo_depth(a_srcfile)
        elif self.r_srffile is not None and self.r_srffile != "":
            hypo_dep = fault_utils.get_hypocenter(a_srffile, sta_base)[2]
        else:
            # No way to determine hypocenter depth, existing
            print("No way to determine hypocenter depth, exiting!")
            sys.exit(1)
        sBBP_in.write("%s\n" % hypo_dep)
        #source model
        sBBP_in.write("1\n")
        #output format 2=velocity
        sBBP_in.write("2\n")
        sBBP_in.flush()
        sBBP_in.close()

        vMname_in = open("%s/VMname.list" % a_tmpdir_mod, "w")
        for stat in site_list:
            vMname_in.write("%s\n" % self.r_velocity)
        vMname_in.flush()
        vMname_in.close()

        #stitch expects <stat>.3comp for LF and <stat>.<comp>.nl1D.001 for HF
        for stat in site_list:
            if self.acc:
                if os.path.exists("%s/%d.%s-lf.acc.bbp" % (a_tmpdir, sim_id, stat.scode)):
                    #integrate to velocity
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s-lf.acc.000 ewfile=%s/%d.%s-lf.acc.090 udfile=%s/%d.%s-lf.acc.ver wcc2bbp=0 < %s/%d.%s-lf.acc.bbp >> %s 2>&1" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, self.log)
                    bband_utils.runprog(cmd)
                    for comp in ['000', '090', 'ver']:
                        cmd = "%s/integ_diff integ=1 filein=%s/%d.%s-lf.acc.%s fileout=%s/%d.%s-lf.vel.%s" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, comp, a_tmpdir, sim_id, stat.scode, comp)
                        bband_utils.runprog(cmd)
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s-lf.vel.000 ewfile=%s/%d.%s-lf.vel.090 udfile=%s/%d.%s-lf.vel.ver units=cm/s wcc2bbp=1 > %s/%d.%s-lf.bbp" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode)
                    bband_utils.runprog(cmd)
                else:
                    print("Can't find LF acceleration file for "
                          "site %s, aborting." % (stat.scode))
                    sys.exit(2)
                if os.path.exists("%s/%d.%s.acc.bbp" % (a_tmpdir, sim_id, stat.scode)):
                    #split into components and rename
                    #nl1D requires 1 line with nt dt
                    #followed by 5 entries per line
                    #Assumes UCSB HF, or why would we be running this?
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s.acc.000 ewfile=%s/%d.%s.acc.090 udfile=%s/%d.%s.acc.ver wcc2bbp=0 < %s/%d.%s.acc.bbp >> %s 2>&1" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, self.log)
                    bband_utils.runprog(cmd)
                    for comp in ['000', '090', 'ver']:
                        cmd = "%s/integ_diff integ=1 filein=%s/%d.%s.acc.%s fileout=%s/%d.%s.vel.%s" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, comp, a_tmpdir, sim_id, stat.scode, comp)
                        bband_utils.runprog(cmd)
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s.vel.000 ewfile=%s/%d.%s.vel.090 udfile=%s/%d.%s.vel.ver units=cm/s wcc2bbp=1 > %s/%d.%s.bbp" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode)
                    bband_utils.runprog(cmd)
                elif os.path.exists("%s/%d.%s-hf.acc.bbp" % (a_tmpdir, sim_id, stat.scode)):
                    print("Using HF only.")
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s-hf.acc.000 ewfile=%s/%d.%s-hf.acc.090 udfile=%s/%d.%s-hf.acc.ver wcc2bbp=0 < %s/%d.%s-hf.acc.bbp >> %s 2>&1" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, self.log)
                    bband_utils.runprog(cmd)
                    for comp in ['000', '090', 'ver']:
                        cmd = "%s/integ_diff integ=1 filein=%s/%d.%s-hf.acc.%s fileout=%s/%d.%s-hf.vel.%s" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, comp, a_tmpdir, sim_id, stat.scode, comp)
                        bband_utils.runprog(cmd)
                    cmd = "%s/wcc2bbp nsfile=%s/%d.%s-hf.vel.000 ewfile=%s/%d.%s-hf.vel.090 udfile=%s/%d.%s-hf.vel.ver units=cm/s wcc2bbp=1 > %s/%d.%s-hf.bbp" % (install.A_GP_BIN_DIR, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode, a_tmpdir, sim_id, stat.scode)
                    bband_utils.runprog(cmd)
                else:
                    print("Can't find HF acceleration file for "
                          "site %s, aborting." % (stat.scode))
            if os.path.exists("%s/%d.%s-lf.bbp" % (a_tmpdir, sim_id, stat.scode)):
                shutil.copy2("%s/%d.%s-lf.bbp" % (a_tmpdir, sim_id, stat.scode), "%s/%d.%s-lf.prestitch.bbp" % (a_tmpdir, sim_id, stat.scode))
                shutil.copy2("%s/%d.%s-lf.bbp" % (a_tmpdir, sim_id, stat.scode), "%s/%s.3comp" % (a_tmpdir, stat.scode))
            else:
                print("Can't find LF velocity file for site %s, aborting." %
                      (stat.scode))
                sys.exit(2)
            if os.path.exists("%s/%d.%s.bbp" % (a_tmpdir, sim_id, stat.scode)):
                #split into components and rename
                #nl1D requires 1 line with nt dt
                #followed by 5 entries per line
                #Assumes UCSB HF, or why would we be running this?
                shutil.copy2("%s/%d.%s.bbp" % (a_tmpdir, sim_id, stat.scode), "%s/%d.%s.prestitch.bbp" % (a_tmpdir, sim_id, stat.scode))
                self.split_bbp(a_tmpdir, stat.scode, "%s/%d.%s.bbp" % (a_tmpdir, sim_id, stat.scode))
            elif os.path.exists("%s/%d.%s-hf.bbp" % (a_tmpdir, sim_id, stat.scode)):
                shutil.copy2("%s/%d.%s-hf.bbp" % (a_tmpdir, sim_id, stat.scode), "%s/%d.%s-hf.prestitch.bbp" % (a_tmpdir, sim_id, stat.scode))
                self.split_bbp(a_tmpdir, stat.scode, "%s/%d.%s-hf.bbp" % (a_tmpdir, sim_id, stat.scode))
            else:
                print("Can't find HF velocity file for site %s, aborting." %
                      (stat.scode))

        cmd = "%s/stitchBBP >> %s 2>&1" % (install.A_UCSB_BIN_DIR, self.log)
        bband_utils.runprog(cmd)

        #reconstitute BBP file
        for stat in site_list:
            print("%s/%s.000.gmBB.001" % (a_tmpdir, stat.scode))
            if os.path.exists("%s/%s.000.gmBB.001" % (a_tmpdir_mod, stat.scode)):
                self.make_bbp(a_tmpdir_mod, stat.scode)
            shutil.copy2("%s/%d.%s-stitch.bbp" %
                         (a_tmpdir_mod, sim_id, stat.scode),
                         "%s/%d.%s.vel.bbp" %
                         (a_tmpdir, sim_id, stat.scode))
            shutil.copy2("%s/%d.%s-stitch.bbp" %
                         (a_tmpdir_mod, sim_id, stat.scode),
                         "%s/%d.%s.vel.bbp" %
                         (a_outdir, sim_id, stat.scode))

            # Create acceleration seismogram

            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, stat.scode))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, stat.scode))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, stat.scode))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.vel.bbp" % (sim_id, stat.scode))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "wcc2bbp=0 < %s >> %s 2>&1" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True)

            for comp in ['000', '090', 'ver']:
                # Differentiate each component
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.%s" %
                                      (sim_id, stat.scode, comp))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.acc.%s" %
                                       (sim_id, stat.scode, comp))

                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/integ_diff diff=1 filein=%s fileout=%s" %
                       (install.A_GP_BIN_DIR, filein, fileout))
                bband_utils.runprog(cmd, abort_on_error=True)

            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.000" % (sim_id, stat.scode))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.090" % (sim_id, stat.scode))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.ver" % (sim_id, stat.scode))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.acc.bbp" % (sim_id, stat.scode))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True)

            # Copy acceleration bbp file to outdir
            shutil.copy2(os.path.join(a_tmpdir, "%d.%s.acc.bbp" %
                                      (sim_id, stat.scode)),
                         os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                      (sim_id, stat.scode)))

        os.chdir(old_cwd)

        print("UCSB Stitch Completed".center(80, '-'))

if __name__ == '__main__':
    VELOCITY = sys.argv[1]
    SRFFILE = sys.argv[2]
    SRCFILE = sys.argv[3]
    STATIONS = sys.argv[4]
    ACC = False
    SIM_ID = 0
    if len(sys.argv) == 7:
        if sys.argv[5] == "True":
            ACC = True
        SIM_ID = int(sys.argv[6])
    else:
        SIM_ID = int(sys.argv[5])
    UC_OBJ = UCStitch(VELOCITY, SRCFILE, SRFFILE, STATIONS,
                      acc=ACC, sim_id=SIM_ID)
    UC_OBJ.run()
