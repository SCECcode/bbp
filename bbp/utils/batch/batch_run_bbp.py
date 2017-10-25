#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# batch_run_bbp.py
# v1.2, 20110819
Process a folder with XML files with broadband platform in batch mode
"""
import optparse
import os
import sys
import shutil
import time
import datetime

from install_cfg import InstallCfg
import bband_utils
import seqnum

# LOG run commands only
LOG_ONLY = True

def findXML():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--in-dir", dest="inDir",
                      help="Input folder with XML files", metavar="IN_DIR")
    parser.add_option("-r", "--resume", action="store_true", dest="resume",
                      help="Resume batch processing of XML files",
                      metavar="RESUME")
    parser.add_option("-o", "--out-dir", dest="outDir",
                      help="Output folder for simulation data (indata, outdata, log, tmpdata) with simID",
                      metavar="OUT_DIR")
    parser.add_option("-f", "--force", action="store_true", dest="force",
                      help="Force overwrite of BBP folders (indata, outdata, log, tmpdata) with same simID",
                      metavar="FORCE")

    (options, args) = parser.parse_args()

    if not options.inDir:
        parser.error("Folder with XML files is required!")

    if os.path.exists(options.inDir):
        files = os.listdir(options.inDir)
    else:
        print "Invalid input dir: %s" % options.inDir
        sys.exit(1)

    if not files:
        print "No XML files were found in input dir: %s" % options.inDir
        sys.exit(1)

    if options.outDir:
        if not os.path.exists(options.outDir):
            try:
                os.mkdir(options.outDir)
            except:
                print ("Failed to create output dir: %s! Aborting..." %
                       (options.outDir))
                sys.exit(1)
        if not os.path.exists("%s/indata" % options.outDir):
            try:
                os.mkdir("%s/indata" % options.outDir)
            except:
                print ("Failed to create dir: %s/indata! Aborting..." %
                       (options.outDir))
                sys.exit(1)
        if not os.path.exists("%s/outdata" % options.outDir):
            try:
                os.mkdir("%s/outdata" % options.outDir)
            except:
                print ("Failed to create dir: %s/outdata! Aborting..." %
                       (options.outDir))
                sys.exit(1)
        if not os.path.exists("%s/tmpdata" % options.outDir):
            try:
                os.mkdir("%s/tmpdata" % options.outDir)
            except:
                print ("Failed to create dir: %s/tmpdata! Aborting..." %
                       (options.outDir))
                sys.exit(1)
        if not os.path.exists("%s/logs" % options.outDir):
            try:
                os.mkdir("%s/logs" % options.outDir)
            except:
                print ("Failed to create dir: %s/logs! Aborting..." %
                       (options.outDir))
                sys.exit(1)

    if options.force:
        print "WARNING: Force overwrite is ON!"
        print "Some existing BBP data folders will be overwritten!"

    install = InstallCfg()
    num_sims = 0
    total_sims = 0
    if options.inDir[-1:] == os.path.sep:
        options.inDir = options.inDir[0:-1]
    files = sorted([f for f in os.listdir(options.inDir) if os.path.isfile(options.inDir + os.path.sep + f)])
    resume_list = ""
    if options.resume==True:
        if os.path.exists("%s/batch_resume.txt" % install.A_OUT_LOG_DIR):
            resume_fp = open("%s/batch_resume.txt" % install.A_OUT_LOG_DIR,
                             'r')
            resume_list = resume_fp.readlines()
            resume_fp.close()
    else:
        if os.path.exists("%s/batch_resume.txt" % install.A_OUT_LOG_DIR):
            os.remove("%s/batch_resume.txt" % install.A_OUT_LOG_DIR)

    run_list = []

    for file in files:
        if file.endswith(".xml"):
            total_sims += 1
            if options.resume == True:
                match = False
                for i in resume_list:
                    if file == i.strip():
                        match = True
                        print "Skipping %s" % file
                        break
                if match == True:
                    continue

            run_list.append(os.path.abspath(os.path.join(options.inDir, file)))
            num_sims += 1

    if not num_sims == total_sims:
        sims = "%d/%d" % (num_sims, total_sims)
    else:
        sims = str(num_sims)

    # Setup the simlist and movelist for logging
    simlist = []
    mvlist = []

    print "Preparing to run %s simulations from %s" % (sims, options.inDir)
    run_count = 1
    for file in run_list:
        filename = os.path.basename(file)
        file_base = filename[0:filename.find(".xml")]
        pieces = file_base.split("_")

        simID =- 1
        if len(pieces) > 1:
            try:
                simID = int(pieces[0])
            except ValueError:
                print "Failed to fetch simID from XML file name: %s" % file

        if simID==-1:
            simID = int(seqnum.getSeqNum())
            print "Running with generated simID: %d" % simID

        t0 = time.time()
        start_time = time.strftime("%Y/%m/%d-%H:%M:%S", time.localtime())

        indatadir = "%s/%d" % (install.A_IN_DATA_DIR, simID)
        outdatadir = "%s/%d" % (install.A_OUT_DATA_DIR, simID)
        tmpdatadir = "%s/%d" % (install.A_TMP_DATA_DIR, simID)
        logdir = "%s/%d" % (install.A_OUT_LOG_DIR, simID)
        logfiledir = "%s/logs/%d" % (options.outDir, simID)
        log_file = "%s/%s.log" % (logfiledir, file_base)
        # Make sure we have an absolute path for log_file
        log_file = os.path.abspath(log_file)
        if not os.path.exists(logfiledir):
            try:
                os.mkdir(logfiledir)
            except:
                print ("Failed to create dir: %s! Aborting..." %
                       (logfiledir))
                sys.exit(1)

        dir_exists = False
        if os.path.exists(indatadir):
            if options.force:
                shutil.rmtree(indatadir)
            else:
                dir_exists = True
        if os.path.exists(tmpdatadir):
            if options.force:
                shutil.rmtree(tmpdatadir)
            else:
                dir_exists = True
        if os.path.exists(outdatadir):
            if options.force:
                shutil.rmtree(outdatadir)
            else:
                dir_exists = True
        if os.path.exists(logdir):
            if options.force:
                shutil.rmtree(logdir)
            else:
                dir_exists = True

        if dir_exists:
            print "BBP folders with simID %d exists!"
            print "Force overwrite is OFF, skipping %s" % (simID, file)
            continue

        print "Processing file (%d/%d): %s" % (run_count, num_sims, file)
        cmd = "%s/run_bbp.py -x %s -s %d -l %s" % (install.A_COMP_DIR,
                                                   file, simID, log_file)
        if (LOG_ONLY == True):
            simlist.append("%s\n" % (cmd))
            if options.outDir:
                # Notes:
                # 1) force option not currently supported while
                # logging sims
                # 2) assumption is that dir outdir/simid
                # does not already exist
                od_indatadir = "%s/indata" % (options.outDir)
                od_outdatadir = "%s/outdata" % (options.outDir)
                od_tmpdatadir = "%s/tmpdata" % (options.outDir)
                od_logdir = "%s/logs" % (options.outDir)
                mvlist.append("mv %s %s\n" % (indatadir, od_indatadir))
                mvlist.append("mv %s %s\n" % (tmpdatadir, od_tmpdatadir))
                mvlist.append("mv %s %s\n" % (outdatadir, od_outdatadir))
                mvlist.append("mv %s %s\n" % (logdir, od_logdir))
            run_count += 1
            continue
        bband_utils.runprog(cmd, False)
        if options.outDir:
            od_indatadir = "%s/indata/%d" % (options.outDir, simID)
            od_outdatadir = "%s/outdata/%d" % (options.outDir, simID)
            od_tmpdatadir = "%s/tmpdata/%d" % (options.outDir, simID)
            od_logdir = "%s/logs/%d" % (options.outDir, simID)

            od_dir_exists = False
            if os.path.exists(od_indatadir):
                if options.force:
                    shutil.rmtree(od_indatadir)
                else:
                    od_dir_exists = True
            if os.path.exists(od_tmpdatadir):
                if options.force:
                    shutil.rmtree(od_tmpdatadir)
                else:
                    od_dir_exists = True
            if os.path.exists(od_outdatadir):
                if options.force:
                    shutil.rmtree(od_outdatadir)
                else:
                    od_dir_exists = True
            if os.path.exists(od_logdir):
                if options.force:
                    shutil.rmtree(od_logdir)
                else:
                    od_dir_exists = True

            if dir_exists:
                print "Warning: Folder(s) with simID %d exists in output folder %s! Force overwrite is OFF, output will be left in BBP folders!" % (simID, options.outDir)
            else:
                if os.path.exists(indatadir):
                    shutil.move(indatadir, od_indatadir)
                if os.path.exists(tmpdatadir):
                    shutil.move(tmpdatadir, od_tmpdatadir)
                if os.path.exists(outdatadir):
                    shutil.move(outdatadir, od_outdatadir)
                if os.path.exists(logdir):
                    shutil.move(logdir, od_logdir)

        run_time = time.time() - t0
        run_time_str = str(datetime.timedelta(seconds=run_time))
        end_time = time.strftime("%Y/%m/%d-%H:%M:%S", time.localtime())

        if options.outDir:
            if not od_dir_exists:
                out_data_dir = od_outdatadir
        else:
            out_data_dir = "%s/%d" % (install.A_OUT_DATA_DIR, simID)

        files = os.listdir(out_data_dir)
        if files:
            try:
                run_log_fp = open("%s/batch_run.log" % options.inDir, 'a')
            except IOError:
                run_log_fp = open("%s/batch_run.log" % install.A_OUT_LOG_DIR,
                                  'a')

            run_log_fp.write("%d\t%s\t%s\t%s\t%s\t%s\n" %
                             (simID, os.path.abspath(file),
                              os.path.abspath(out_data_dir),
                              start_time, end_time, run_time_str))
            run_log_fp.flush()
            run_log_fp.close()
            resume_fp = open("%s/batch_resume.txt" % install.A_OUT_LOG_DIR,
                             'a')
            resume_fp.write("%s\n" % filename)
            resume_fp.flush()
            resume_fp.close()
        run_count +=1

    # Write the sims to the execlog
    print "Opening %s/batch_run_bbp_sims.log" % (options.inDir)
    execlog = open("%s/batch_run_bbp_sims.log" % (options.inDir), 'w')
    for sim in simlist:
        execlog.write(sim)
    execlog.close()

    # Write the moves to the execlog
    print "Opening %s/batch_run_bbp_moves.log" % (options.inDir)
    execlog = open("%s/batch_run_bbp_moves.log" % (options.inDir), 'w')
    for mv in mvlist:
        execlog.write(mv)
    execlog.close()

if __name__ == '__main__':
    findXML()
