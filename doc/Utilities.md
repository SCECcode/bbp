# Broadband Platform Utilities

This page contains a brief description of the utilities included in the Broadband Platform distribution.

### Plotting velocity models

The 'plot_velocity_model.py' tool included in the 'utils/misc' directory can be used to create 1D velocity profile plots of the velocity models included with the Broadband Platform. It uses as input a velocity model profile in the format specified in the [File Formats](./File-Format-Guide). In addition to the output plot file (in PNG format), users can also specify a title for the plot. For example:

```
 $ plot_velocity_model.py -v nr02-vs500.fk1d -o labasin500.png -t "LA Basin 500"
```

### Multi-segment scripts

The Broadband Platform includes 5 methods that are able to run simulations containing a rupture composed of multiple planes or segments. The GP, SDSU, SONG, Irikura Recipe Method 1, and Irikura Recipe Method 2 methods are able to read multiple simple source description (SRC) files, each containing information about one rupture segment, and create a multi-segment, version 2.0 SRF file that can be used by other modules in the Broadband Platform to compute synthetic seismograms. When running a multi-segment validation event with these methods, users do not need to do anything different than what they do for single segment events.

It is also possible to calculate ground motions from a multi-segment rupture one segment at a time. This can be used to calculate the contribution of each segment to the overall ground motions, or during the development and debugging of simulation methods. In order to do this, users will need to run each segment separately and then merge the results of all segments into a multi-segment simulation. This is done with two merge scripts available in 'utils/misc', merge_multisegment_validation.py and merge_multisegment_scenario.py.

Users can use the '--help' flag to see the command-line options in more detail but, basically, users will need to provide a list of the simulation IDs of each of the individual segments so the script will read the data and combine the results. For example, in the 1992 Landers earthquake model with 3 separate segments, users will first need to run the Broadband Platform 3 times (once for each of the segments). With multi-segment validation packages the Platform will ask the user which segment the user wants to use. Then, results can be merged with the following command:

```
 $ merge_multisegment_validation.py --sim-id 104 --event LandersMS 101 102 103
```

In the example above, we assume that Landers segments 1, 2, and 3 were calculated and their simulation ids are 101, 102, and 103, respectively. The merge tool will combine the simulations under simulation id 104.

### Cluster scripts

The 'bbp_discovery_validation.py' and 'bbp_discovery_scenario.py' tools, available in 'utils/batch' are used to run Broadband Platform validation and scenario simulations on the USC's CARC cluster. The scripts are specific to the USC cluster and use the Slurm job scheduler and file systems specific to it. It is our hope that users can adapt these scripts without too much effort to work on other computing environments.

Each script includes several command-line options that allow users to customize their simulations but, basically, for validation simulations, users will need to specify a simulation method, the validation event, the number of realizations (or simulation variations) to run, whether to use the site response module or not, a simulation directory, certain randomization parameters, and an e-mail address so the user is notified when the job starts and finishes. For example:

```
 $ bbp_discovery_validation.py -s -c gp -e tottori -n 64 --email <user_email> --no-hypo-rand -d gp-tottori-64r
```

The above command creates 64 realizations for the Tottori event using the 'GP' simulation method. It includes the site response module in the simulations and disables the randomization of the hypocenter (each of the 50 simulations will still include a different slip distribution and random seed). The tool will create 50 realizations of the simulation in the 'gp-tottori-64r' folder and will print the Slurm command users need to type to submit the simulation to the cluster.

For scenario simulations, instead of providing a validation event users provide a source description file, a station list, and a simulation region. For example:

```
 $ bbp_discovery_scenario.py -c exsim -n 64 -v LABasin500 --src <source_desc.src> --stl <station_list.stl> --email <user_email> --hypo-rand -d exsim-my-rup-64r
```

Similarly, the tool will generate 64 realizations for the source description and station list specified in the command-line flags in the 'exsim-my-rup-64r' and will print the Slurm command users need to type to submit the simulation to the cluster.

Each script includes several options that allow users to customize the simulation, please use the '--help' flag for a complete list of the available options.

### Comparing realizations

When used with a cluster validation simulation, the 'compare_bbp_runs.py' tool (available in 'utils/batch') can be used to compare the different realizations and rank them on how well each one matches the observed data. It takes only the top-level simulation directory (as specified with the '-d' in the cluster scripts above). For example:

```
 $ compare_bbp_runs.py gp-nr-64r
 Simulation    Average Bias
  10000022    0.077551
  10000006    0.099381
  10000028    0.102009
  10000002    0.102098
  ...
  ...
  ...
  10000015    0.185119
  10000033    0.185664
  10000020    0.192202
```

In the output above, the 64 realizations in 'gp-nr-64r' are ranked from best to worst based on their Average Bias. In this case, the best realization is number 10000022 with a 0.077551 average bias and the worst is realization 10000020 with an average bias of 0.192202.

### Multi-segment cluster runs

Multi-segment simulations outside of the cluster environment were described above. For cluster simulations, the same validation and scenario cluster scripts used for single-segment cluster simulations can also be used to run multi-segment runs. Multi-segment simulations include multiple source description (SRC) files, one for each rupture segment. At this time on the Broadband Platform, there are two methodologies for calculating seismograms from multi-segment ruptures.

#### Computing ground motions for all segments together

The GP, SDSU, SONG, Irikura Recipe Method1, and Irikura Recipe Method 2 methods compute a single version 2.0 extended rupture SRF file containing data from all SRC files. When simulations complete, they contain the final result of all segments combined, and nothing else is needed. For these two methods, please use the same cluster scripts listed above, without any changes to the usage.

#### Calculating segments separately

For the methods above, it is also possible to compute the segments separately. In this case, each rupture segment needs to be calculated separately and all segments need to be merged later. For example, the multi-segment version of the 1992 Landers earthquake contains 3 segments. They can be calculated with the following commands:

```
 $ bbp_discovery_validation.py -c gp -e "Landers Multi Segment" -n 64 --email <user_email> --no-hypo-rand -d gp-landers1-64r -s --segment 1
```

```
 $ bbp_discovery_validation.py -c gp -e "Landers Multi Segment" -n 64 --email <user_email> --no-hypo-rand -d gp-landers2-64r -s --segment 2 --variation 2 --first-seg gp-landers1-64r
```

```
 $ bbp_discovery_validation.py -c gp -e "Landers Multi Segment" -n 64 --email <user_email> --no-hypo-rand -d gp-landers3-64r -s --segment 3 --variation 3 --first-seg gp-landers1-64r
```

With each of the command above, the Broadband Platform will create 64 realizations for each of the 3 rupture segments in the 1992 Landers earthquake. Users will be prompted to type the 'sbtach' command to submit each of these simulations to the cluster. Please note that simulations for segments 2 and 3 use information in the segment 1 simulation (the random seed in the first segments' source description file) so that certain randomized simulation parameters (e.g. rupture velocity) remain consistent across all segments.

After all simulation complete, users will need to merge the results for the three segments into a single Broadband Platform simulation, using the 'bbp_merge_multisegment_validation.py' tool available in the 'utils/batch' folder. For the Landers example above with the 3 segments, users can type:

```
 $ bbp_merge_multisegment_validation.py -d gp-landers-ms-64r -e "Landers Multi Segment" gp-landers1-64r gp-landers2-64r gp-landers3-64r
```

For each of the 64 realizations, the command above will merge the synthetic seismograms and then compute the post-processing analysis steps. This process can take several hours.

### Convergence

The 'bbp_converge.py' tool included in the 'utils/batch' directory can be used with the other cluster utilities to generate a summary showing how many realizations are required for a certain method to converge. This can be useful for estimating how many simulations (realizations) need to be used for each specific simulation method. The 'input_directory' is the same top-level directory used by the cluster scripts to generate the simulations (specified with the '-d' flag), and the 'method' parameter should be the simulation method used for the simulations (e.g. gp, sdsu, etc) and is used only for the plot title.

```
 $ bbp_converge.py --ns 1000 -i <input_directory> -o <output_plot_file.png> -c <method>
```

### Importing Simulation Data

The "import_bbp_simulation.py" tool included in the "utils/misc" folder can be used to import a set of simulation results calculated outside of the BBP into the Platform for post-processing. For example, users can compute simulation results using a method outside of the BBP and then run the standard BBP post-processing tools to compare their results against recorded data from one of the BBP validation events. For a list of validation events available on the Platform, please refer to folders inside the $BBP_VAL directory.

The first step in running the Broadband Platform to process data calculated outside of the BBP is to import the simulated data and convert it to the format expected by the BBP. Users can use the "import_bbp_simulation.py" tool available in the "utils/misc" directory to perform some of these steps. The tool assumes the seismograms are already in BBP format but need to be renamed and copied into a simulation directory.

The command-line options available for the "import_bbp_simulation.py" tool are:

```
 $ import_bbp_simulation.py --help
 usage: import_bbp_simulation.py [-h] --sim_id SIM_ID --src_file SRC_FILE
                                --station_list STATION_LIST --input_dir
                                INPUT_DIR [--prefix PREFIX] [--suffix SUFFIX]

 Create a BBP simulation from a set of seismograms.

 optional arguments:
  -h, --help            show this help message and exit
  --sim_id SIM_ID       simulation id
  --src_file SRC_FILE, --src SRC_FILE
                        src file
  --station_list STATION_LIST, --stl STATION_LIST
                        station list
  --input_dir INPUT_DIR, -i INPUT_DIR
                        input directory
  --prefix PREFIX, -p PREFIX
                        prefix for input files
  --suffix SUFFIX, -s SUFFIX
                        suffix for input files
```

The tool assumes users have a folder containing velocity and acceleration files for a set of stations. It provides some flexibility on the naming but velocity and acceleration filenames will need to follow the format:

```
 <prefix><station_name><suffix>.vel.bbp and <prefix><station_name><suffix>.acc.bbp
```

The prefix and suffix can be specified using command-line options and are optional. The station_names need to match the station names from the station list. The tool will parse the station list and look for both velocity and acceleration files for each station -- both are required. If users have simulated data for all stations for a validation event, they can simply use the station list file already available in the validation package. If data is available for only a subset of stations, users will first need to make a copy of the station list and edit it to delete stations not included in the simulation set.

Once users have set up their station list (or confirmed that the station list from the validation event matches the data that is available), users can run the "import_bbp_simulation.py" tool. The tool will copy the input files and create a simulation folder for this simulation.

As an example, assuming the Landers validation event, an input folder named "my_sims" and a "p-" prefix and a "-model1" suffix, users can run the tool with:

```
 $ import_bbp_simulation.py --prefix "p-" --suffix="-model1" --input_dir my_sims
   --src_file $BBP_VAL/Landers/common/landers_v14_02_1.src
   --station_list $BBP_VAL/Landers/common/landers_v19_06_2.stl --sim_id 20200501
```

Please note the equal and extra quotes for the suffix parameter as the name starts with a dash. Also, please make sure the BBP_DATA_DIR variable is set up before running the tool so that the files can be copied to the right place. The value provided for the sim_id needs to be unique (not one from a previously ran simulation).

Once the simulation files are copied, users will need to run the BBP. This will be a two-step process. First, users will need to create an XML file for the simulation. Users can run the BBP with the following command-line options:

```
 $ run_bbp.py -g --expert
```

The "-g" flag tells the BBP not to run the simulation itself, only generate the XML file. The "--expert" flag may be used to customize the simulation parameters (e.g. selecting a custom station file). Here's one example of running the Broadband Platform with the Landers validation event:

```
 $ run_bbp.py --expert -g
 Welcome to the SCEC Broadband Platform version 22.4.0.
 ================================================================================

 Please select the Broadband Platform mode of operation:
   * Validation - Simulates a historical event
   * Scenario   - Runs a user-defined hypothetical event

 Do you want to perform a validation simulation (y/n)? y
 ================================================================================

 Please select a validation event from the list below:

 (1) Alum Rock
 (2) Chino Hills
 (3) Landers
 (4) NR
 (5) North Palm Springs
 (6) Parkfield
 (7) SanSimeon
 (8) Whittier
 ? 3
 ================================================================================

 The Broadband Platform includes several scientific methods that can be used
 to calculate synthetic seismograms.

 Choose a Method to use in this Broadband validation simulation:
 (1) GP (Graves & Pitarka)
 (2) UCSB
 (3) SDSU
 (4) EXSIM
 (5) Song
 (6) Irikura Recipe Method 1 (Irikura1)
 (7) Irikura Recipe Method 2 (Irikura2)
 ? 1
 ================================================================================

 Each validation package includes a default source description (SRC) file for a
 historical event. Would you like to provide a different file instead of the default
 file provided? Answer 'no' here if you would like to use the standard source
 file for this event.

 Do you want to provide a custom source file (y/n)? n
 ================================================================================
 SRC file: /home/jane/bbp/bbp_val/Landers/common/landers_v14_02_1.src
 ================================================================================

 When starting a simulation from a source description (SRC) file, the Broadband
 Platform workflow should include a rupture generator. Answer 'yes' here unless
 providing a complex Standard Rupture Format (SRF) file.

 Do you want to run the rupture generator (y/n)? y
 ================================================================================

 Station Selection
 =================
 Would you like to:
   (1) generate seismograms for all stations in the validation package
       OR
   (2) provide a custom list with a subset of the stations
 ? 1
 ================================================================================
 STL file: /home/jane/bbp/bbp_val/Landers/common/landers_v19_06_2.stl
 ================================================================================

 Site Response
 =============
 Running a site response module is an optional step while running a Broadband
 Platform simulation. It requires a station list file containing the Vs30 values
 for each station location.

 Do you want to run the site response module (y/n)? y

 Please select a site response module to use in this Broadband simulation:
 (1) GP (Graves & Pitarka)
 (2) PySeismoSoil
 ? 1
 ================================================================================
 Do you want to generate velocity seismograms' plots (y/n)? y
 ================================================================================
 Do you want to generate acceleration seismograms' plots (y/n)? y
 ================================================================================

 The Broadband Platform can generate comparison plots of the validation data
 against GMPEs to show how GMPEs match the recorded data for a certain event.

 Do you want to generate a GMPE comparison plot (y/n)? y
 ================================================================================

 Please select a GMPE set to use in the comparison (number of name are ok):

 (1) CENA GROUP 1
 (2) NGA-West1
 (3) NGA-West2
 ? 3
 ================================================================================

 Goodness-of-Fit Plot
 ====================
 Running a goodness-of-fit (GoF) module is an optional  step while running a
 Broadband Platform simulation. It creates a comparison plot showing how well the
 calculated seismograms fit recorded data.

 Do you want to run a goodness-of-fit module (y/n)? y
 ================================================================================

 Users can optionally select a Goodness of Fit module to plot a comparison of
 how well the simulated seismograms match the recorded data in a historical event.

 Choose a Goodness of Fit (GOF) Module:
 (1) GP
 (2) FAS
 (3) GP and SDSU
 (4) GP and FAS
 ? 1
 ================================================================================

 Additional Metrics
 ==================
 Calculating additional metrics is an optional step on the Broadband Platform.
 It creates additional plots and data files that can be used to study the simulation.

 Do you want to calculate additional metrics (y/n)? n
 XML file /home/jane/bbp/bbp_data/xml/9996390.xml generated, not running workflow.
```

Please note the xml file generated by the BBP (it is listed in the last line). A couple of things to note are that since the simulation data is already generated outside of the Platform, it does not matter what simulation method is selected (not what site response module is selected). If users created a custom station list, they should select the "custom" option after selecting to run the rupture generator and then provide the path for the station list created for this simulation. Finally, users should select the GP Goodness-of-Fit module and request velocity and acceleration plots.

Once that step is completed, users will need to run the BBP a second time using the following options:

```
 $ run_bbp.py -m -s <sim_id> -x <xml_file> -r plot_map
```

The "sim_id" should match the value specified by the user when running the "import_bbp_simulation.py" tool. The "xml_file" should be the xml file created by the BBP when it was called the first time (see above). It is important not to forget to specify the "-r plot_map" flag as this flag will tell the BBP to skip the seismogram generation steps and go directly to the post-processing stage.
