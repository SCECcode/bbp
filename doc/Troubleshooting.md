If you experience trouble building the platform or successfully running test and simulations, try the following solutions.

### Build Errors

The instruction for installing Broadband Platform are listed in Section 1: "Installing the Second-Generation Broadband Platform" of the User guide. If after following all steps listed in the Broadband Platform Installation page the build fails, check if the failure is listed in this section and try the solution to fix the issue you are facing.

### Unit/Acceptance Test Warnings and Errors

Unit and Acceptance tests are provided to verify the Broadband platform and it's supporting modules built by the user are functioning as designed. Under certain circumstances, some of these test might fail or print a warning message. While some of these failures might indicate serious problems that will have to be addressed before the platform can used, it is acceptable to ignore some of the warnings. This section lists some issues and their solutions.

#### test_bbtoolbox

The test_bbtoolbox unit test will fail on Mac OS X. This is due to the raytracing code in the SDSU method not working on Mac OS X. This is a known issue and we hope to have this issue resolved in a future version of the Broadband Platform.

#### test_gof

In certain environments, the test_gof unit test will display a series of warning messages regarding the availability of certain fonts used in the plots. The user would see a number of messages like:

```
 /usr/lib64/python2.7/site-packages/matplotlib/font_manager.py:1242: UserWarning: findfont: Font family ['STIXSizeThreeSym'] not found.
 Falling back to Bitstream Vera Sans
 /usr/lib64/python2.7/site-packages/matplotlib/font_manager.py:1242: UserWarning: findfont: Font family ['STIXSizeFourSym'] not found.
 Falling back to Bitstream Vera Sans
```

These messages do not indicate a failure of the unit test. At the end of this test you will still see the "ok" status, indicating that the test was successful.

#### test_user-SDSU
#### test_valid-northridge-SDSU
#### test_valid-northridge-SDSU_seis

All acceptance tests for the SDSU method will fail on Mac OS X due to the raytracing code in the SDSU method. This is a known issue and we hope to have this resolved in an upcoming version of the Broadband Platform.

### Runtime Warnings and Errors

Below is a list of known issues in the current Broadband 16.5.0 release that can sometimes cause a simulation to fail.

#### SDSU Method: BBToolbox exits with code 174

The current version of SDSU's BBToolbox code can sometimes fail with an exit code of 174. This will cause the Broadband simulation to stop immediately and a message similar to the one below will be recorded in the bbtoolbox log file (which is located in the simulation's log directory at $BBP_DATA_DIR/logs/<simulation_id>:

```
 *** Welcome to the Broad-Band Toolbox (v1.6) ***
 
 Initialising the code, please type input filename ...
 rl=          131
 rl=          131
 Starting slug3d: by J. Vidale, 1988, UCSC
 Starting punch: by J. Hole, 1993, UBC-Stanford
 -- travel-times written to :time3d_P.out
 -- selected velocity file  :vel3d_P.bin
 WARNING: Computing only to max radius = 1001
 Completed radius = 10
 Completed radius = 20
 Completed radius = 30
 Completed radius = 40
 Completed radius = 50
 Completed radius = 60
 Completed radius = 70
 Completed radius = 80
 Completed radius = 90
 Completed radius = 100
 Completed radius = 110
 wavefront done 
 Starting slug3d: by J. Vidale, 1988, UCSC
 Starting punch: by J. Hole, 1993, UBC-Stanford
 -- travel-times written to :time3d_S.out
 -- selected velocity file  :vel3d_S.bin
 WARNING: Computing only to max radius = 1001
 forrtl: severe (174): SIGSEGV, segmentation fault occurred
 Image              PC                Routine            Line        Source             
 BBtoolbox.exe      000000000043A28B  Unknown               Unknown  Unknown
 BBtoolbox.exe      0000000000403A8B  Unknown               Unknown  Unknown
 BBtoolbox.exe      000000000040372C  Unknown               Unknown  Unknown
 libc.so.6          000000337641ECDD  Unknown               Unknown  Unknown
 BBtoolbox.exe      0000000000403629  Unknown               Unknown  Unknown
```

This message reflects an issue in the raytracer code used by BBToolbox that will cause the code to fail with certain geometries. The SDSU team is currently working on a fix to this issue and we expect a future version of the Platform to include the fix. Meanwhile, the user can try to change the problem geometry in order to avoid the error. This can be done by overriding the geometry that is calculated based on the stations and rupture (using a larger x/y/z grid will not affect the simulation results but can sometimes solve this raytracer issue). To see what is the current problem geometry, users can look for the following line in the .bbpar file (located in the $BBP_DATA_DIR/indata/<simulation_id> directory:

```
 /* GRID DEFINITION [X-Y-Z] FOR RAYTRACING: "FAR-SIDE" (IN KM) */
 130.0 140.0 125.0
```

To specify a different geometry to be used by BBToolbox's raytracer, users can include the following keys in the SRC file:

```
 GRID_X = xxx
 GRID_Y = yyy
 GRID_Z = zzz
```

Where xxx, yyy, and zzz are the new dimensions for the raytracer to use (specified in km).

#### SDSU Method: BBToolbox exits with code 139

The current version of SDSU's BBToolbox code can sometimes fail with an exit code of 139. We are working with the modelers to resolve this issue, but it appears related to memory management issues in the code. We have found that it is sometimes useful to increase the stack size (if one is set). This should be done before the simulation starts, with the command below (for bash):

```
 $ ulimit -s unlimited
```

#### SDSU Method: BBToolbox fails unit test with *** buffer overflow detected *** error

The current version of SDSU's BBToolbox code has an internal issue that in some systems trigger built-in compiler checks and cause the code to abort. This is known to happen on more recent Ubuntu/XUbuntu installations, but can also happen in other systems where these checks are enabled by default. When running the BBToolbox unit test users will see a message like:

```
 *** buffer overflow detected ***: /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin/BBtoolbox.exe terminated
 ======= Backtrace: =========
 /lib/x86_64-linux-gnu/libc.so.6(+0x7338f)[0x2afe8fe5138f]
 /lib/x86_64-linux-gnu/libc.so.6(__fortify_fail+0x5c)[0x2afe8fee8c9c]
 /lib/x86_64-linux-gnu/libc.so.6(+0x109b60)[0x2afe8fee7b60]
 /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin/BBtoolbox.exe[0x44351b]
 /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin/BBtoolbox.exe[0x4020e4]
 /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin/BBtoolbox.exe[0x4036bd]
 /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf5)[0x2afe8fdffec5]
 /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin/BBtoolbox.exe[0x401d79]
```

Until this problem is fixed, one solution is disabling the built-in check. That can be done by editing the makefile in the /home/sarah/bbp/16.5.0/bbp/src/sdsu/bin directory so that the first line looks like:

```
 UFLAGS = -w -Wall -U_FORTIFY_SOURCE
```

and the line under the BBtoolbox.exe line looks like:

```
 ${FC} ${MODULES} ${CODES} -O3 -U_FORTIFY_SOURCE ray3DJHfor.o -o BBtoolbox.exe
```

Then, type:

```
 $ make clean
 $ make
```

#### UCSB Method: Site Response module

In the current version of the Broadband Platform, the UCSB site response module is producing unusually high amplitudes. This is likely due to a unit mismatch between the new versions of the source and wave propagation codes, and the old version of the site response module. Therefore, the UCSB method site response module should not be used in this release. We are currently working on a fix and will have a solution implemented by the next release of the Platform.