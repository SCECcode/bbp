# Broadband Platform Utilities

This page contains a brief description of the utilities included in the Broadband Platform distribution.

### Plotting velocity models

The 'plot_velocity_model.py' tool included in the 'utils/misc' directory can be used to create 1D velocity profile plots of the velocity models included with the Broadband Platform. It uses as input a velocity model profile in the format specified in the [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide). In addition to the output plot file (in PNG format), users can also specify a title for the plot. For example:

```
 $ plot_velocity_model.py -v nr02-vs500.fk1d -o labasin500.png -t "LA Basin 500"
```

### Multi-segment scripts

The Broadband Platform includes 4 methods that are able to run simulations containing a rupture composed of multiple planes or segments. The SONG and Irikura Recipe Method 1 methods are able to read multiple simple source description (SRC) files, each containing information about one rupture segment, and create a multi-segment, version 2.0 SRF file that can be used by other modules in the Broadband Platform to compute synthetic seismograms. When running a multi-segment validation event with these two methods, users do not need to do anything different than what they do for single segment events.

The GP and SDSU methods are also able to calculate ground motions for multi-segment ruptures. However, at this time, this is done one segment at a time. Users will need to run each segment separately and then merge the results of all segments into a multi-segment simulation. This is done with two merge scripts available in 'utils/misc', merge_multisegment_validation.py and merge_multisegment_scenario.py.

Users can use the '--help' flag to see the command-line options in more detail but, basically, users will need to provide a list of the simulation IDs of each of the individual segments so the script will read the data and combine the results. For example, in the 1992 Landers earthquake model with 3 separate segments, users will first need to run the Broadband Platform 3 times (once for each of the segments). With multi-segment validation packages the Platform will ask the user which segment the user wants to use. Then, results can be merged with the following command:

```
 $ merge_multisegment_validation.py --sim-id 104 --event LandersMS 101 102 103
```

In the example above, we assume that Landers segments 1, 2, and 3 were calculated and their simulation ids are 101, 102, and 103, respectively. The merge tool will combine the simulations under simulation id 104.

### Cluter scripts

The 'bbp_hpcc_validation.py' and 'bbp_hpcc_scenario.py' tools, available in 'utils/batch' are used to run Broadband Platform validation and scenario simulations on the USC's HPC cluster. The scripts are specific to the USC cluster and use the Slurm job scheduler and file systems specific to it. It is our hope that users can adapt these scripts without too much effort to work on other computing environments.

Each script includes several command-line options that allow users to customize their simulations but, basically, for validation simulations, users will need to specify a simulation method, the validation event, the number of realizations (or simulation variations) to run, whether to use the site response module or not, a simulation directory, certain randomization parameters, and an e-mail address so the user is notified when the job starts and finishes. For example:

```
 $ bbp_hpcc_validation.py -s -c gp -e tottori -n 50 --email <user_email> --no-hypo-rand -d gp-tottori-50r
```

The above command creates 50 realizations for the Tottori event using the 'GP' simulation method. It includes the site response module in the simulations and disables the randomization of the hypocenter (each of the 50 simulations will still include a different slip distribution and random seed). The tool will create 50 realizations of the simulation in the 'gp-tottori-50r' folder and will print the Slurm command users need to type to submit the simulation to the cluster.

For scenario simulations, instead of providing a validation event users provide a source description file, a station list, and a simulation region. For example:

```
 $ bbp_hpcc_scenario.py -c exsim -n 50 -v LABasin500 --src <source_desc.src> --stl <station_list.stl> --email <user_email> --hypo-rand -d exsim-my-rup-50r
```

Similarly, the tool will generate 50 realizations for the source description and station list specified in the command-line flags in the 'exsim-my-rup-50r' and will print the Slurm command users need to type to submit the simulation to the cluster.

Each script includes several options that allow users to customize the simulation, please use the '--help' flag for a complete list of the available options.

### Comparing realizations

When used with a cluster validation simulation, the 'compare_bbp_runs.py' tool (available in 'utils/batch') can be used to compare the different realizations and rank them on how well each one matches the observed data. It takes only the top-level simulation directory (as specified with the '-d' in the cluster scripts above). For example:

```
 $ compare_bbp_runs.py gp-nr-50r
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

In the output above, the 50 realizations in 'gp-nr-50r' are ranked from best to worst based on their Average Bias. In this case, the best realization is number 10000022 with a 0.077551 average bias and the worst is realization 10000020 with an average bias of 0.192202.

### Multi-segment cluster runs

Multi-segment simulations outside of the cluster environment were described above. For cluster simulations, the same validation and scenario cluster scripts used for single-segment cluster simulations can also be used to run multi-segment runs. Multi-segment simulations include multiple source description (SRC) files, one for each rupture segment. At this time on the Broadband Platform, there are two methodologies for calculating seismograms from multi-segment ruptures.

#### SONG and Irikura Recipe Method 1

The SONG and Irikura Recipe Method 1 methods compute a single version 2.0 extended rupture SRF file containing data from all SRC files. When simulations complete, they contain the final result of all segments combined, and nothing else is needed. For these two methods, please use the same cluster scripts listed above, without any changes to the usage.

#### Other simulation methods

Other simulation methods that support multi-segment simulations are the GP and SDSU methods. For these two methods however, each rupture segment needs to be calculated separately and all segments need to be merged later. For example, the multi-segment version of the 1992 Landers earthquake contains 3 segments. They can be calculated with the following commands:

```
 $ bbp_hpcc_validation.py -c gp -e "Landers Multi Segment" -n 50 --email <user_email> --no-hypo-rand -d gp-landers1-50r -s --segment 1
```

```
 $ bbp_hpcc_validation.py -c gp -e "Landers Multi Segment" -n 50 --email <user_email> --no-hypo-rand -d gp-landers2-50r -s --segment 2 --variation 2 --first-seg gp-landers1-50r
```

```
 $ bbp_hpcc_validation.py -c gp -e "Landers Multi Segment" -n 50 --email <user_email> --no-hypo-rand -d gp-landers3-50r -s --segment 3 --variation 3 --first-seg gp-landers1-50r
```

With each of the command above, the Broadband Platform will create 50 realizations for each of the 3 rupture segments in the 1992 Landers earthquake. Users will be prompted to type the 'sbtach' command to submit each of these simulations to the cluster. Please note that simulations for segments 2 and 3 use information in the segment 1 simulation (the random seed in the first segments' source description file) so that certain randomized simulation parameters (e.g. rupture velocity) remain consistent across all segments.

After all simulation complete, users will need to merge the results for the three segments into a single Broadband Platform simulation, using the 'bbp_merge_multisegment_validation.py' tool available in the 'utils/batch' folder. For the Landers example above with the 3 segments, users can type:

```
 $ bbp_merge_multisegment_validation.py -d gp-landers-ms-50r -e "Landers Multi Segment" gp-landers1-50r gp-landers2-50r gp-landers3-50r
```

For each of the 50 realizations, the command above will merge the synthetic seismograms and then compute the post-processing analysis steps. This process can take several hours.

### Convergence

The 'bbp_converge.py' tool included in the 'utils/batch' directory can be used with the other cluster utilities to generate a summary showing how many realizations are required for a certain method to converge. This can be useful for estimating how many simulations (realizations) need to be used for each specific simulation method. The 'input_directory' is the same top-level directory used by the cluster scripts to generate the simulations (specified with the '-d' flag), and the 'method' parameter should be the simulation method used for the simulations (e.g. gp, sdsu, etc) and is used only for the plot title.

```
 $ bbp_converge.py --ns 1000 -i <input_directory> -o <output_plot_file.png> -c <method>
```
