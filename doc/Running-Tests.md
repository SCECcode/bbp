The Broadband Platform contains two kinds of tests: unit and acceptance tests.

Unit tests are designed to test each module separately, using a set of input files, and compare the results against known outputs. They verify that each module has been built and is working correctly. Acceptance tests verify that the modules are working correctly together as a simulation method. They test the Platform end-to-end using different combinations with known inputs and compare the results. Each method (for both user-defined and validation events) are tested as checks against integration errors. Both unit and acceptance tests require the LABasin Greens' Functions package to be installed. The acceptance tests also require users to install the Northridge (NR) validation package.

### Running Unit Tests

In order to run the unit tests, users should go to the tests directory and type:

```
 $ cd /home/sarah/bbp/19.4.0/bbp/tests
```

We recommend you run the Unit and Acceptance tests in the background using a command like:

```
 nohup ./UnitTests.py > $BBP_DATA_DIR/tests.log 2>&1 &
```

This will output stdout to a file called tests.log. The UnitTests also produce a log file in the BBP_DATA_DIR with additional details about each method tested. Examine the tests.log file to determine whether the tests complete successfully. When they do you should see a message like the following:

You can run the Unit tests directly, with the following command. However you must keep the terminal window open until the tests complete.

```
 $ ./UnitTests.py
```

The tests should begin and will take between 20-60 minutes to run, depending on your computer speed. At the end of each test, a "ok" should be printed if the test was successful. At the end, the program will print the number of tests that passed and the number of tests that failed. If a test has failed, first check that you have built the executables. You can rerun just the specific test that failed (test_<module>.py). If the test is still failing, also verify that you have the ref_data directory, since it contains the input and reference files. If you're looking for more information about the failure, you can consult the Unit Tests log file in `$BBP_DATA_DIR/logs/unit_tests.log`. If you can't determine the reason for the failure, contact software [at] scec.org for support.

Please note that one of the SDSU components, BBToolbox, does not pass its Unit Test when running on a Mac OS X system and the Platform will automatically skip this test. This is a known limitation of the current BBP release. You will see a message similar to:

```
 *** Mac OS X detected: skipping SDSU BBToolbox unit test.
```

This reduces the available number of ground motion methods that can be run on the Mac. On Linux, there are seven methods available (GP, SDSU, UCSB, EXSIM, SONG, Irikura Recipe Method 1, Irikura Recipe Method 2). On a Mac, there are six methods available (GP, UCSB, EXSIM, SONG, Irikura Recipe Method 1, Irikura Recipe Method 2).

Once the unit tests all pass, proceed to the acceptance tests. If there are any failure or errors while running the unit tests, consult the Troubleshooting section at the end of this user guide for know issues and their solutions. Unit test results on a Mac look like this:

```
 $ nohup ./UnitTests.py > $BBP_DATA_DIR/tests.log 2>&1 &
 $ more $BBP_DATA_DIR/tests.log
 test_runprog (test_bband_utils.TestBBandUtils) ... ok
 test_runprog2 (test_bband_utils.TestBBandUtils) ... ok
 test_runprog3 (test_bband_utils.TestBBandUtils) ... ok
 test_execute_platform_bbp (test_python_code.TestPythonCode) ... ok
 test_python_code_comps (test_python_code.TestPythonCode) ... ok
 test_python_code_tests (test_python_code.TestPythonCode) ... ok
 test_arias_duration (test_arias.TestArias) ... ok
 test_bbp2peer (test_bbp_format.TestBBPFormat) ... ok
 test_exsim2bbp (test_bbp_format.TestBBPFormat) ... ok
 test_peer2bbp (test_bbp_format.TestBBPFormat) ... ok
 test_vm2ucsb (test_vm2vm.TestVm2vm) ... ok
 test_vm2ucsb_nga (test_vm2vm.TestVm2vm) ... ok
 test_vm2vm (test_vm2vm.TestVm2vm) ... ok
 test_vm2vm_nga (test_vm2vm.TestVm2vm) ... ok
 test_xy2ll (test_cc.TestCC) ... ok
 test_genslip (test_genslip.TestGenslip) ... ok
 test_jbsim (test_jbsim.TestJbsim) ... ok
 test_hfsims (test_hfsims.TestHfsims) ... ok
 test_wcc_siteamp (test_wcc_siteamp.TestWccSiteamp) ... ok
 test_match (test_match.TestMatch) ... ok
 test_gensrf (test_gensrf.TestGenSRF) ... ok
 test_irikura_hf (test_irikura_hf.TestIrikuraHF) ... ok
 test_uc_fault_utils (test_uc_fault_utils.TestUCFaultUtils) ... ok
 test_ucgen (test_ucrmg.TestUCrmg) ... ok
 test_syn1d (test_syn1d.TestSyn1D) ... ok
 test_uc_site (test_uc_site.TestUCSite) ... ok
 test_bbtoolbox (test_bbtoolbox.TestBBToolbox) ... ok
 test_exsim (test_exsim.TestExsim) ... ok
 test_rmg (test_rmg.TestRMG) ... ok
 test_rotd50 (test_rotd50.TestRotD50) ... ok
 test_rotd100 (test_rotd100.TestRotD100) ... ok
 test_gof (test_gp_gof.TestGPGof) ... ok
 test_mogof_stat (test_sdsu_mogof.TestSDSUMOGof) ... ok
 test_anderson_gof (test_anderson_gof.TestAndersonGof) ... ok
 test_rzz2015 (test_rzz2015.TestRZZ2015) ... ok
 test_as16 (test_as16.TestAS16) ... ok

----------------------------------------------------------------------
Ran 36 tests in 1455.138s

OK
```

### Running Acceptance Tests

Make sure the unit tests pass before moving on to the acceptance tests. To run the acceptance tests, users should:

```
 $ cd /home/sarah/bbp/19.4.0/bbp/tests
```

To run the acceptance tests in a terminal background process, on a Linux computer, the following command will run the tests as a background process, and redirect stderr and stdout to a "accept_tests.log" file.

```
 $ nohup ./AcceptTests.py > $BBP_DATA_DIR/accept_tests.log 2>&1 &
```

You can run the Acceptance tests directly in a terminal window, but you must keep the terminal window open until all tests complete.

```
 $ ./AcceptTests.py
```

There will be one test for each method in the platform, for both the validation mode (historical events) and the scenario modes (hypothetical earthquakes). These tests will take somewhere between 6 and 24 hours to run, depending on how fast the computer is.

When they're complete, the console will either print "ok" or how many tests failed. Acceptance test failures indicate that the modules are not integrated correctly. Like with the unit tests, verify that you have the ref_data directory.  If a certain acceptance test fails, you can get more information by consulting the acceptance test logs in `$BBP_DATA_DIR/logs/acceptance_test_logs/<test that failed>.log`. If you can't determine the reason for the failure, contact software [at] scec.org for support.

Here is an example output from the tests.log file created by running this command line above.

```
-bash-4.1$ more $BBP_DATA_DIR/accept_tests.log
 nohup: ignoring
 ==> Number of tests to run: 14
 test_user-EXSIM (__main__.BBPAcceptanceTests) ... ok
 test_user-GP (__main__.BBPAcceptanceTests) ... ok
 test_user-IRIKURA_RECIPE_M1 (__main__.BBPAcceptanceTests) ... ok
 test_user-IRIKURA_RECIPE_M2 (__main__.BBPAcceptanceTests) ... ok
 test_user-SDSU (__main__.BBPAcceptanceTests) ... ok
 test_user-SONG (__main__.BBPAcceptanceTests) ... ok
 test_user-UCSB (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-EXSIM (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-GP (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-IRIKURA_RECIPE_M1 (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-IRIKURA_RECIPE_M2 (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-SDSU (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-SONG (__main__.BBPAcceptanceTests) ... ok
 test_valid-northridge-UCSB (__main__.BBPAcceptanceTests) ... ok

 ----------------------------------------------------------------------
 Ran 14 tests in 25981.848s

 OK
```

Please note that since the SDSU method will produce errors on Mac OS X systems, two of the tests above will be skipped. This is expected on a Mac OS X computer.
