#! /usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Tests for the full validations in the Broadband Platform.
$Id: test_full_validations.py 1734 2016-09-13 17:38:17Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import optparse
import unittest

# Import Broadband modules
import seqnum
import bband_utils
from install_cfg import InstallCfg

# Global variables
SIM_ID = None

class ValidationTest(unittest.TestCase):
    """
    Validation test case for the Broadband Platform
    """
    def validation_test(self, event, code_base, test_name, run_rupture_gen):
        """
        This function runs the full validation test for the event and
        code_base specified. It creates a Broadband option file with
        the required answers in the input directory and then invokes
        the platform. After execution finishes, it checks if the .bbp
        files were correctly generated, and looks for the gof output
        file.
        """
        # event - Northridge = 1, Loma Prieta = 2, Landers = 3, LOMAP = 4
        #         NR = 5, Whittier = 6, Saguenay = 7
        # code_base GP = 1, UCSB = 2, SDSU = 3
        events = ["Northridge", "Loma Prieta",
                  "Landers", "LOMAP", "NR", "Whittier", "Saguenay"]

        num_stations = 0
        bbp_files_found = 0
        gof_png_file_found = False
        gof_png_file_valid = False
        a_station_list_filename = None
        install = InstallCfg()

        if SIM_ID is None:
            self.sim_id = int(seqnum.get_seq_num())
        else:
            self.sim_id = SIM_ID

        #print
        #print "Using sim_ID: %d..." % self.sim_id
        # Create input and log directories
        cmd = "mkdir -p %s" % os.path.join(install.A_IN_DATA_DIR,
                                           str(self.sim_id))
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Error creating input directory")
        cmd = "mkdir -p %s" % os.path.join(install.A_OUT_LOG_DIR,
                                           str(self.sim_id))
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Error creating log directory")
        # Create path for option file
        a_opt_filename = os.path.join(install.A_IN_DATA_DIR,
                                      str(self.sim_id),
                                      "%d.optfile" % (self.sim_id))
        try:
            opt_file = open(a_opt_filename, 'w')
            # Select validation run
            opt_file.write("y\n")
            # Select event
            opt_file.write("%s\n" % (events[event - 1]))
            # Select the rupture generator
            if run_rupture_gen:
                # Yes to rupture generator
                opt_file.write("y\n")
                if code_base == 2:
                    # Run UCSB rupture generator
                    opt_file.write("2\n")
                else:
                    # Otherwise run GP rupture generator
                    opt_file.write("1\n")
                # No need to provide custom src file, use existing
                opt_file.write("n\n")
            else:
                # Not using the rupture generator
                opt_file.write("n\n")
            # Validate all stations
            opt_file.write("1\n")
            # Select codebase
            if code_base == 1:
                # GP across the board
                opt_file.write("1\n1\n4\n")
            elif code_base == 2:
                # UCSB low-frequency requires UCSB across the board
                opt_file.write("3\n2\n4\n")
            elif code_base == 3:
                # SDSU doesn't have a low-frequency component
                opt_file.write("1\n3\n4\n")
            else:
                self.failIf(True, "Code base not supported!")
            # We don't want velocity seismograms
            opt_file.write("y\n")
            # But we want acceleration seismograms
            opt_file.write("y\n")
            # Select GOF module, for validation we want to use GP
            opt_file.write("1\n")
            opt_file.close()
        except IOError:
            self.failIf(True, "Cannot create option file!")

        # Log output to a_log_filename
        a_log_filename = os.path.join(install.A_OUT_LOG_DIR, str(self.sim_id),
                                      "%s.%s.log" % (str(self.sim_id),
                                                     test_name))

        # Run BBP command
        cmd = "%s/run_bbp.py -s %d -o %s >> %s" % (install.A_COMP_DIR,
                                                   self.sim_id,
                                                   a_opt_filename,
                                                   a_log_filename)
        # print cmd
        # Run the validation scenario
        self.failIf(bband_utils.runprog(cmd, False) != 0,
                    "Fatal error while running the Broadband Platform")
        # First look for the station list
        try:
            files = os.listdir(os.path.join(install.A_IN_DATA_DIR,
                                            str(self.sim_id)))
        except OSError:
            self.failIf(True, "Cannot access input directory!")

        for my_file in files:
            if my_file.endswith(".stl"):
                a_station_list_filename = os.path.join(install.A_IN_DATA_DIR,
                                                       str(self.sim_id),
                                                       my_file)
                break

        # Fail if we don't have a station list
        self.failIf(a_station_list_filename is None,
                    "Cannot find file with station list!")

        # Now, figure out how many stations we have
        try:
            stl_file = open(a_station_list_filename, 'r')
            for _ in stl_file:
                num_stations = num_stations + 1
        except IOError:
            self.failIf(True, "Cannot read station list file!")

       # Now look for the output files
        try:
            files = os.listdir(os.path.join(install.A_OUT_DATA_DIR,
                                            str(self.sim_id)))
        except OSError:
            self.failIf(True, "Cannot access output directory!")
        for my_file in files:
            if my_file.endswith(".acc.bbp"):
                a_bbp_file_path = os.path.join(install.A_OUT_DATA_DIR,
                                               str(self.sim_id),
                                               my_file)
                try:
                    my_size = os.stat(a_bbp_file_path).st_size
                except OSError:
                    self.failIf(True,
                                "Cannot stat file %s!" % a_bbp_file_path)
                # Ok, let's count this file
                if my_size > 0:
                    bbp_files_found = bbp_files_found + 1
            if my_file.startswith("gof-") and my_file.endswith(".png"):
                # Found the gof.png file
                gof_png_file_found = True
                a_gofpng_file_path = os.path.join(install.A_OUT_DATA_DIR,
                                                  str(self.sim_id),
                                                  my_file)
                try:
                    my_size = os.stat(a_gofpng_file_path).st_size
                except OSError:
                    self.failIf(True,
                                "Cannot stat file %s!" % a_gofpng_file_path)
                # Check if size > 0
                if my_size > 0:
                    gof_png_file_valid = True

        # Now check if we got the right results
        self.failIf(not gof_png_file_found, "GOF file not found!")
        self.failIf(not gof_png_file_valid, "GOF file size is 0!")
        self.failIf(bbp_files_found != num_stations,
                    "Found %d valid BBP files, expecting %d!" %
                    (bbp_files_found, num_stations))

    def test_northridge_gp(self):
        """
        Validation test for the Northridge GP scenario
        """
        self.validation_test(event=1, code_base=1,
                             test_name="test_northridge_gp",
                             run_rupture_gen=False) # Northridge, GP

    def test_northridge_ucsb(self):
        """
        Validation test for the Northridge UCSB scenario
        """
        self.validation_test(event=1, code_base=2,
                             test_name="test_northridge_ucsb",
                             run_rupture_gen=False) # Northridge, UCSB

    def test_northridge_sdsu(self):
        """
        Validation test for the Northridge SDSU scenario
        """
        self.validation_test(event=1, code_base=3,
                             test_name="test_northridge_sdsu",
                             run_rupture_gen=False) # Northridge, SDSU

    def test_lp_gp(self):
        """
        Validation test for the Loma Prieta GP scenario
        """
        self.validation_test(event=2, code_base=1,
                             test_name="test_lp_gp",
                             run_rupture_gen=False) # Loma Prieta, GP

    def test_lp_ucsb(self):
        """
        Validation test for the Loma Prieta UCSB scenario
        """
        self.validation_test(event=2, code_base=2,
                             test_name="test_lp_ucsb",
                             run_rupture_gen=False) # Loma Prieta, UCSB

    def test_lp_sdsu(self):
        """
        Validation test for the Loma Prieta SDSU scenario
        """
        self.validation_test(event=2, code_base=3,
                             test_name="test_lp_sdsu",
                             run_rupture_gen=False) # Loma Prieta, SDSU

    def test_la_gp(self):
        """
        Validation test for the Landers GP scenario
        """
        self.validation_test(event=3, code_base=1,
                             test_name="test_la_gp",
                             run_rupture_gen=False) # Landers, GP

    def test_la_ucsb(self):
        """
        Validation test for the Landers UCSB scenario
        """
        self.validation_test(event=3, code_base=2,
                             test_name="test_la_ucsb",
                             run_rupture_gen=False) # Landers, UCSB

    def test_la_sdsu(self):
        """
        Validation test for the Landers SDSU scenario
        """
        self.validation_test(event=3, code_base=3,
                             test_name="test_la_sdsu",
                             run_rupture_gen=False) # Landers, SDSU

    def test_lomap_gp(self):
        """
        Validation test for the LOMAP GP scenario
        """
        self.validation_test(event=4, code_base=1,
                             test_name="test_lomap_gp",
                             run_rupture_gen=True) # LOMAP, GP

    def test_lomap_ucsb(self):
        """
        Validation test for the LOMAP UCSB scenario
        """
        self.validation_test(event=4, code_base=2,
                             test_name="test_lomap_ucsb",
                             run_rupture_gen=True) # LOMAP, UCSB

    def test_lomap_sdsu(self):
        """
        Validation test for the LOMAP SDSU scenario
        """
        self.validation_test(event=4, code_base=3,
                             test_name="test_lomap_sdsu",
                             run_rupture_gen=True) # LOMAP, SDSU

    def test_nr_gp(self):
        """
        Validation test for the NR GP scenario
        """
        self.validation_test(event=5, code_base=1,
                             test_name="test_nr_gp",
                             run_rupture_gen=True) # NR, GP

    def test_nr_ucsb(self):
        """
        Validation test for the NR UCSB scenario
        """
        self.validation_test(event=5, code_base=2,
                             test_name="test_nr_ucsb",
                             run_rupture_gen=True) # NR, UCSB

    def test_nr_sdsu(self):
        """
        Validation test for the NR SDSU scenario
        """
        self.validation_test(event=5, code_base=3,
                             test_name="test_nr_sdsu",
                             run_rupture_gen=True) # NR, SDSU

    def test_whittier_gp(self):
        """
        Validation test for the WHITTIER GP scenario
        """
        self.validation_test(event=6, code_base=1,
                             test_name="test_whittier_gp",
                             run_rupture_gen=True) # WHITTIER, GP

    def test_whittier_ucsb(self):
        """
        Validation test for the WHITTIER UCSB scenario
        """
        self.validation_test(event=6, code_base=2,
                             test_name="test_whittier_ucsb",
                             run_rupture_gen=True) # WHITTIER, UCSB

    def test_whittier_sdsu(self):
        """
        Validation test for the WHITTIER SDSU scenario
        """
        self.validation_test(event=6, code_base=3,
                             test_name="test_whittier_sdsu",
                             run_rupture_gen=True) # WHITTIER, SDSU

    def test_saguenay_gp(self):
        """
        Validation test for the SAGUENAY GP scenario
        """
        self.validation_test(event=7, code_base=1,
                             test_name="test_saguenay_gp",
                             run_rupture_gen=True) # SAGUENAY, GP

    def test_saguenay_ucsb(self):
        """
        Validation test for the SAGUENAY UCSB scenario
        """
        self.validation_test(event=7, code_base=2,
                             test_name="test_saguenay_ucsb",
                             run_rupture_gen=True) # SAGUENAY, UCSB

    def test_saguenay_sdsu(self):
        """
        Validation test for the SAGUENAY SDSU scenario
        """
        self.validation_test(event=7, code_base=3,
                             test_name="test_saguenay_sdsu",
                             run_rupture_gen=True) # SAGUENAY, SDSU
# -----------------------------------------------------------------------------
# Main starts here
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # Parse command-line option
    prog_base = os.path.split(sys.argv[0])[1]
    prog_usage = ("usage: %s [options]" % prog_base)
    prog_desc = ("User needs to specify CODE_BASE and EVENT_ID,"
                 " or use the --all option")
    parser = optparse.OptionParser(usage=prog_usage, description=prog_desc)
    parser.add_option("-s", "--simID", dest="input_sim_id", type="int",
                      help="Force a simID", metavar="SIM_ID")
    parser.add_option("-a", "--all", action="store_true", dest="run_all",
                      help="Runs all validation tests, cannot be used with -s")
    parser.add_option("-e", "--event", action="store", type="string",
                      dest="event", metavar="EVENT_ID",
                      help="Event, valid options are: northridge, lp, la, "
                           "lomap, nr, whittier, saguenay")
    parser.add_option("-c", "--code", action="store", type="string",
                      dest="codebase", metavar="CODE_BASE",
                      help="Codebase to use, options are: GP, SDSU, and UCSB")
    (options, args) = parser.parse_args()


    # Check options...
    if options.input_sim_id and options.run_all:
        parser.error("Options -a and -s cannot be used at the same time.")

    if options.input_sim_id:
        try:
            SIM_ID = int(options.input_sim_id)
        except ValueError:
            parser.error("SIM_ID must be a positive integer")
        if SIM_ID < 1:
            parser.error("SIM_ID must be a positive integer")
    else:
        SIM_ID = None

    if (options.run_all is None and
        (options.event is None or options.codebase is None)):
        # Must pick codebase and event unless select "all"
        parser.print_help()
        sys.exit(1)

    # Check events
    if options.event is not None:
        options.event = options.event.lower()
        if (options.event != "northridge" and
            options.event != "lp" and
            options.event != "lomap" and
            options.event != "la" and
            options.event != "nr" and
            options.event != "whittier" and
            options.event != "saguenay"):
            parser.error("Event must be: northridge, lp, la, lomap, nr, whittier, or saguenay")

    # Check codebase
    if options.codebase is not None:
        options.codebase = options.codebase.lower()
        if (options.codebase != "gp" and
            options.codebase != "ucsb" and
            options.codebase != "sdsu"):
            parser.error("Codebase must be: gp, ucsb, or sdsu")

    # Instantiate test suite
    suite = unittest.TestSuite()

    if options.run_all:
        # If --all, let's create all tests
        suite.addTest(ValidationTest('test_northridge_gp'))
        suite.addTest(ValidationTest('test_northridge_ucsb'))
        suite.addTest(ValidationTest('test_northridge_sdsu'))
        suite.addTest(ValidationTest('test_lp_gp'))
        suite.addTest(ValidationTest('test_lp_ucsb'))
        suite.addTest(ValidationTest('test_lp_sdsu'))
        suite.addTest(ValidationTest('test_la_gp'))
        suite.addTest(ValidationTest('test_la_ucsb'))
        suite.addTest(ValidationTest('test_la_sdsu'))
        suite.addTest(ValidationTest('test_lomap_gp'))
        suite.addTest(ValidationTest('test_lomap_ucsb'))
        suite.addTest(ValidationTest('test_lomap_sdsu'))
        suite.addTest(ValidationTest('test_nr_gp'))
        suite.addTest(ValidationTest('test_nr_ucsb'))
        suite.addTest(ValidationTest('test_nr_sdsu'))
        suite.addTest(ValidationTest('test_whittier_gp'))
        suite.addTest(ValidationTest('test_whittier_ucsb'))
        suite.addTest(ValidationTest('test_whittier_sdsu'))
        suite.addTest(ValidationTest('test_saguenay_gp'))
        suite.addTest(ValidationTest('test_saguenay_ucsb'))
        suite.addTest(ValidationTest('test_saguenay_sdsu'))
    else:
        # Otherwise, we just pick the test the user selected
        if options.event == "northridge":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_northridge_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_northridge_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_northridge_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "lp":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_lp_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_lp_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_lp_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "la":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_la_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_la_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_la_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "lomap":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_lomap_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_lomap_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_lomap_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "nr":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_nr_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_nr_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_nr_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "whittier":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_whittier_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_whittier_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_whittier_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        elif options.event == "saguenay":
            if options.codebase == "gp":
                suite.addTest(ValidationTest('test_saguenay_gp'))
            elif options.codebase == "ucsb":
                suite.addTest(ValidationTest('test_saguenay_ucsb'))
            elif options.codebase == "sdsu":
                suite.addTest(ValidationTest('test_saguenay_sdsu'))
            else:
                # This should not happen!
                parser.error("Invalid code option!")
        else:
            # This should not happen!
            parser.error("Invalid code option!")

    # Now run test(s)
    unittest.TextTestRunner(verbosity=2).run(suite)
