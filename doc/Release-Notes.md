### Broadband Platform 19.4.0

This full release of the Broadband Platform includes the following features and bug fixes. Below is a summary of the improvements and modifications includes in this release of the Broadband Platform.

#### Method Updates

* Added the Irikura Recipe Method 2 simulation method to the Broadband Platform. It uses the same rupture generator as the Irikura Recipe Method 1, developed by Arben Pitarka, followed by the GP low-frequency wave propagation code. It then uses a high-frequency module developed by NIED.

* Both SONG and Irikura Recipe Method 1 methods now support multi-segment ruptures. These two methods create a version 2.0 SRF file that can be used directly by the GP wave propagation codes and allows a single Broadband Platform simulation to include all segments from a multi-plane rupture.

* The GP site response module is now used in all California and Japan validation events.

* Re-calculated GFs for California and Japan using 500m/s as reference Vs, bringing the computed synthetic seismograms closer to site conditions and requiring smaller corrections by the site amplification module.

#### Other Improvements

* Updated the GP rupture generator to genslip-5.4.1 and the GP high-frequency code to hb_high_6.0.3.

* Added a Central California simulation region with its own velocity model, GFs, and simulation parameters for all available methods.

* GP site response module is now available for use with all 7 simulation methods.

* Upgraded cluster simulation scripts for use with the Slurm job scheduler. The original collection of scripts using PBS is still available for reference.

### Broadband Platform 17.3.0

This full release of the Broadband Platform includes the following features and bug fixes. Below is a summary of the improvements and modifications includes in this release of the Broadband Platform.

#### Method Updates

* Added the Irikura Recipe Method 1 simulation method to the Broadband Platform. It includes a new rupture generator module developed by Arben Pitarka to generate a SRF file. Then, it uses the GP method wave propagation codes to generate low- and high-frequency seismograms.

* Updated GP rupture generator code to genslip-5.2.2. This new version of genslip accepts three new parameters used for the simulation of multi-segment ruptures, where time series can be generated segment by segment and then added together before post-processing is done.

* Updated GP match.py module to improve the merging of low- and high-frequency seismograms.

#### General Improvements

* Added a FAS calculation module to the Platform. Work contributed by Jeff Bayless using David Boore's 'smc2fs2' and 'asc2smc' codes. The FAS module produces per-station plots with both N/S and E/W components, along with the smoothed EAS (effective amplitude spectrum) of the two horizontal components.

* Included the calculation of zeta parameter in the RZZ2015 module.

* Added the Central United States simulation region, contributed by Mehrdad Hosseini and Paul Somerville. This includes a new set of Green's functions, calculated up to 1800km.

* Fixed issue in the Anderson GoF and RZZ2015 codes that was causing time series to not align correctly. Thanks to Kim Olsen and Rumi Takedatsu for reporting this bug.

* Added two new scripts: 'merge_multisegment_validation.py' and 'merge_multisegment_scenario.py' to the utils/misc directory. These scripts can be used to combine time series from a number of separately-calculated segments, allowing for a multi-segment rupture to be simulated.

#### Cluster Improvements

* bbp_hpcc_validation.py and bbp_hpcc_scenario.py include option for user to override default walltime. This allows users to specify a short walltime if they know a job will run quickly, allowing the job to be potentially quickly scheduled by PBS.

* Added a '-s' option to both bbp_hpcc_validation.py and bbp_hpcc_scenario.py to enable the use of the site response module. Currently only the GP site response module is supported (it is used for all methods).

* Several modifications to the cluster scripts to enable multi-segment ruptures to be simulated. Added a '--segment' option to specify the segment number, a '--variation' option to enable the use of different sets of random seeds in the cluster. Also added a '--firstsegment' option to enable the scripts to find the first segment of a multi-segment run so that common seeds can be used for certain parameters across multiple segments (used by the GP method).

* Cluster scripts now save metadata file on top-level simulation directory including all command-line options used to generate the cluster simulation. This is useful to track simulation parameters.

* Added two scripts, 'bbp_merge_multisegment_validation.py' and 'bbp_merge_multisegment_scenario.py' to the utils/batch directory. These scripts work similarly to the ones in utils/misc but can be used to calculate multi-segment runs using segments calculated on the cluster.

### Broadband Platform 16.5.0

This full release of the Broadband Platform includes several new features and bug fixes. Here is a summary of the Trac items describing what is new in this release:

* Trac #298: Define evaluation criteria for comparing BBP Releases
* Trac #302: Fix matplotlib deprecation warnings for versions 1.3 and greater
* Trac #303: Integrate RMG rupture generator module into BBP Platform
* Trac #304: Update GP rupture generator code to version 5.0
* Trac #305: Update BBToolbox to version 1.6
* Trac #309: Update GP Rupture Generator to version 5.0.1
* Trac #310: Update GP site response module to version 2014
* Trac #312: Fix compilation warnings for various codes in the GP method
* Trac #313: Modify cluster scripts to use scec-00 filesystem
* Trac #314: Include parameters from RZZ2015 into BBP workflows
* Trac #315: Include parameters from AS2016 in the BBP Workflows
* Trac #318: Add PyProj version to software_status file
* Trac #325: Respect module fails to run for observed seismograms
* Trac #326: Add Anderson GOF to Broadband Platform
* Trac #328: Update LABasin, Mojave and NoCal GFs
* Trac #329: Add RotD100/RotD50 ratio metric to BBP

### Broadband Platform 15.3.0

This full release of the Broadband Platform includes several new features and bug fixes. Here is a summary of the Trac items describing what is new in this release:

* Trac #158 - Platform should write simulation results on BBP_DATA_DIR
* Trac #200 - Add version number to BBP XML workflow description file
* Trac #207 - Use real station distances in platform
* Trac #227 - Convert UCSB method to use GNU compilers
* Trac #228 - Convert SDSU method to use GNU compilers
* Trac #239 - Add the NGA-West2 GMPE models to the Broadband Platform
* Trac #261 - Make simulation inputs available in data output directory
* Trac #262 - Add configurable line width option for BBP seismogram plots
* Trac #263 - Include 3 GMPE models in the Broadband Platform for use in the East
* Trac #264 - Enable LF seismograms module to accept input directory
* Trac #265 - Update Match module to work with user-provided seismograms
* Trac #266 - Update station list generator from PyNGA to allow for footwall and hanging wall stations
* Trac #267 - Remove display requirement from station generation script
* Trac #270 - Fix header in rd50 files in outdata/obs_seis directory in GMPE runs
* Trac #271 - Integrate the UCSB 14.8 code into BBP
* Trac #274 - Add option to validation cluster script to run only rupture generator
* Trac #275 - Allow users to see missing velocity models when selecting a validation event
* Trac #276 - Add random SEED used in simulations to HTML summary in outdata
* Trac #277 - Copy SRC and method-specific parameter files to outdata
* Trac #278 - Integrate SDSU BBToolbox v1.5.5.1 into Broadband Platform
* Trac #279 - Run BBP on cluster using the same command on all nodes
* Trac #280 - Add rule to GP method to calculate DX/DY automatically
* Trac #281 - Label BBP methods that are not yet validated
* Trac #282 - Update SDSU BBtoolbox to fix bug when calculating seismograms at >70km
* Trac #283 - Modify directory structure in validation packages to add common input files directory
* Trac #285 - Error parsing XML file created by BBP
* Trac #286 - ExSim results incorrect for some close by stations
* Trac #287 - Removed unused or obsolete codes and data files from Platform
* Trac #288 - Increase array size in RotD50 code to avoid truncating timeseries
* Trac #289 - Integrate Jan 2015 version of the UCSB source and wave propagation code
* Trac #290 - Update BBP distribution directory structure
* Trac #291 - Improve BBP command-line interface with more intuitive questions
* Trac #292 - Rename LOMAP velocity model to NOCAL
* Trac #293 - Include method name on html report in the outdata directory
* Trac #294 - Update SDSU codebase with latest February 2015 changes
* Trac #295 - Add support for varying rupture velocity in GP method
* Trac #299 - README references icc compilers
* Trac #300 - Add expert mode to run_bbp.py
* Trac #301 - Fix units tests test_vm2vm and test_s2v

### Broadband Platform 14.3.0

This full release of the Broadband Platform includes several new features and bug fixes. Here is a summary of the Trac items describing what is new in this release:

* Trac # 114 - Convert from f77 to gfortran
* Trac # 146 - Move DT out of src file
* Trac # 195 - Produce GMPE plot for validation simulations
* Trac # 212 - wrong unit on respect plot
* Trac # 213 - Remove unneeded GP_GOF parameters
* Trac # 214 - Reduce volume of logged data
* Trac # 215 - Modify SRF plot so that it uses same ALONG_STRIKE/DOWN_DIP reference as SRC file
* Trac # 219 - SRC file should specify Mw, platform then calculates the seismic moment
* Trac # 222 - Integrate new GP Rupture Generator v3.2.1
* Trac # 224 - Add randomization in BBToolbox using the iseed parameter
* Trac # 230 - Modify station list parser to accept float vs30 values
* Trac # 232 - Add the decimation factor patch from SDSU to BBToolbox
* Trac # 233 - Show epicenter location on all map plots
* Trac # 234 - Add '08 suffix to current GMPE plot labels
* Trac # 235 - Modify GoF plot labels
* Trac # 237 - Rotate all multi-figure plots so they are viewed in landscape
* Trac # 240 - Run LF and HF components in the UCSB method using separate GFs
* Trac # 241 - Modify SDSU BBtoolbox to handle smaller-magnitude events
* Trac # 242 - Integrate GP Rupture generator v3.3 into the Broadband Platform
* Trac # 243 - Update GP jbsim LF program
* Trac # 244 - Arias duration module fails when processing zeroed observation files
* Trac # 245 - Make str_fac a region-specific parameter in BBToolbox
* Trac # 246 - Capture and store resource information from shell
* Trac # 247 - Allow user to specify hypocenter randomization area in cluster script
* Trac # 248 - Add option to cluster script for user to run only the rupture generator
* Trac # 249 - Add option to cluster script for users to specify SRF files to use
* Trac # 250 - Correct NGA model labels on GMPE box plots
* Trac # 251 - Integrate ExSIM14 into Broadband Platform
* Trac # 252 - Integrate new version of GP high frequency code
* Trac # 253 - Simplify makefiles and remove support for user-specified compilers
* Trac # 254 - Integrate Feb 2014 version of UCSB code into Broadband Platform
* Trac # 256 - Integrate SDSU model version 1.5.4.1 into the Broadband Platform
* Trac # 257 - Use SEED in SRC file for high frequency simulation in GP method
* Trac # 258 - Require same parameters for cluster post-processing scripts
* Trac # 259 - Produce map plots for GMPE runs
* Trac # 260 - Update SDSU BBToolbox code to version 1.5.5

### Broadband Platform 13.9.0

This full release of the Broadband Platform includes several new features and bug fixes. Here is a summary of the Trac items describing what is new in this release:

* Trac # 38 - Running validation events with SDSU and incorrect station names causes a crash
* Trac # 57 - Setting fmax=20.0 Hz in bbtoolbox_cfg.py generates NaNs in output bbp
* Trac # 80 - Add NGA attenuenation relationship
* Trac # 82 - Separate validation files from GF
* Trac # 83 - Separate GF into individual events
* Trac # 84 - User selectable GF
* Trac # 87 - Integrate Atkinson module into broadband
* Trac # 88 - Integrate Irikura Module into BBP
* Trac # 89 - Provide no site correction option
* Trac # 94 - Add number cutoff distance and number of stations to bias plots
* Trac # 96 - Add link to bias plot results into validation table
* Trac # 98 - Rename main run script to run_bbp.py
* Trac # 100 - Rename urs modules to gp modules
* Trac # 101 - add test to check for mixed spaces and tabs
* Trac # 102 - Move plot files out of Greens Functions
* Trac # 103 - Make site response optional in the platform
* Trac # 104 - Use environment variables to find Greens Functions and Validation directories
* Trac # 106 - Review BBP bbtoobox unit tests on broadband.usc.edu
* Trac # 112 - Collect additional metadata using CSEP environment script and generate full manifest for software distribution
* Trac # 120 - Link velocity model and code to specific gf
* Trac # 121 - full validation runs with src file or user selected srf file
* Trac # 122 - Allow users to provide an alternative directory for input/output data files
* Trac # 123 - reliable system for self reporting software version for broadband
* Trac # 124 - Plot of station map with SDSU code base
* Trac # 125 - Migrate pbs script and parallel scripts into svn trunk
* Trac # 126 - move plots into bbp home directory
* Trac # 130 - Amp Fac Unit test fails due to long filenames
* Trac # 131 - Make parallel scripts configurable to run a number of concurrent instances
* Trac # 136 - Move y2r2b.cpt file from URS_DATA/plot to BBP distribution plot/data directory
* Trac # 137 - Rupture model png not copied to outdata directory
* Trac # 140 - Add Arias Duration plots for bbp full validation sims
* Trac # 142 - Add rotd50 routine to BBP
* Trac # 143 - URS validation fails due to long file paths
* Trac # 144 - Automatically adjust time series plot window to capture entire event
* Trac # 147 - Convert peer obs files to bbp format obs
* Trac # 148 - lowfreq corner of -99 in station list
* Trac # 149 - Check in uwo EXSIM code into trunk with tests
* Trac # 150 - Include the rupture generator in the validation workflow
* Trac # 152 - Rupture plot appears distorted
* Trac # 153 - Generate html option copies extra files to output directory
* Trac # 154 - check rotd50 code
* Trac # 155 - Fix units label on full validation seismogram plots
* Trac # 156 - Enable rotD50 module to handle bbp and peer inputs files
* Trac # 157 - Validate Broadband using ifort version 12.0 on HPCC
* Trac # 159 - Make FMAX user-configurable in BBToolbox
* Trac # 160 - Use RotD50 data in Bias plot instead of average horizontals
* Trac # 161 - Bias plot should only include data within band-pass filter
* Trac # 162 - Show band-pass bars on per-station plots
* Trac # 163 - Create new bias plot showing fit over station distance
* Trac # 164 - Move md5sums into validation and gf packages
* Trac # 165 - Derive UCSB station list from Broadband station list
* Trac # 166 - Derive BBToolbox inputs from Broadband inputs for validation runs
* Trac # 167 - The 2 horizontal components in the PSa5/RotD50 bias plot are inverted
* Trac # 168 - Create a KML file with stations and fault line
* Trac # 169 - Extend RotD50 GOF plot range 0.1Hz - 100Hz
* Trac # 171 - The geobb_srf script fails to find fault corners in certain scenarios
* Trac # 172 - Check station names' length and abort if above max limit
* Trac # 173 - Package GMPE code into the Broadband Platform
* Trac # 174 - Introduce version numbers for validation and GFs packages
* Trac # 175 - Integrate scripts to produce GMPE boxplots into the Platform
* Trac # 176 - Randomize hypocenter location when running multiple validation realizations
* Trac # 177 - Use same UCSB rupture generator binary for vertical and dipping faults
* Trac # 178 - Plot station map with fault using SRC file
* Trac # 179 - Create single component GOF plot
* Trac # 180 - Generate SRF file in XYZ format for the SDSU method
* Trac # 181 - Add Qp and Qs to the SDSU velocity model file
* Trac # 182 - Make switch for randomizing hypocenter required on cluster script
* Trac # 183 - Make SDSU seismograms module take regular station list
* Trac # 185 - Implement resume workflow feature on the Broadband Platform
* Trac # 184 - Check if SRF file exists if user wants to skip rupture generator
* Trac # 186 - Make BBToolbox use magnitude from SRC file
* Trac # 187 - Calculate HYPO_DEPTH automatically for the UCSB method
* Trac # 188 - Add option for users to run the Broadband Platform on the background
* Trac # 189 - Optimize arias_duration script
* Trac # 190 - Remove duplicate tests in UnitTest.py
* Trac # 191 - Add simulation timestamp to index.html file
* Trac # 192 - Integrate UNR Composite Source Model into the Platform
* Trac # 194 - Enable cluster to run user-defined simulations
* Trac # 196 - Create a map GOF plot with color
* Trac # 197 - Create combined GOF for all realizations and stations
* Trac # 198 - Whenever running binary files, print an error message if it is not found
* Trac # 199 - Print error when velocity model doesn't exist for selected method/event
* Trac # 201 - Update the BBToolbox version in the Platform to 1.5
* Trac # 202 - PlotMap.py doesn't work for Matplotlib v1.2.0
* Trac # 203 - Integrate updated version of UCSB rupture generator
* Trac # 204 - Support multiple colorsets on GOF plots
* Trac # 206 - Regenerate acceptance tests for trunk
* Trac # 208 - Add check to make sure number of stations is under UCSB's syn1D limit
* Trac # 209 - Fix race condition for matplotlib cache file when running on the cluster
* Trac # 210 - Increase array size is respect in order to handle larger seismograms
* Trac # 211 - Increase CSM's station limit and check if station list is under new limit

Note: Trac tickets refer to the software bugs/features described in the [BBP Trac Site](https://northridge.usc.edu/trac/broadband)
