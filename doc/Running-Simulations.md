The platform supports two kinds of simulations, validation events and user-defined events. Validation simulations are performed using a historical event, and are directly compared to observed seismograms using goodness-of-fit. User-defined events are run using a rupture description provided by the user which may not necessarily be a historical earthquake.

When you run a simulation, the platform assigns an ID to it. This ID can be used to track the simulation and locate the output data products.

To supply input files to the platform, put them in the run directory. This directory is located in the $BBP_DATA_DIR directory and it is created by the platform when it is first ran. Extensions are important - the platform recognizes station lists (.stl), SRF files (.srf), and simple source descriptions (.src). For example, when selecting a station file, the Broadband Platform will list all the .stl files in the run directory and will prompt the user to select one.

### Platform Modes

The platform can be run in multiple modes. The default is interactive mode, in which the user is prompted to answer a series of questions. Once all the information has been gathered, the simulation run begins.

For a large number of runs, or if the user is repeating a specific run, this can be tedious. The platform provides two other ways to describe a run, with an option file or an XML description.

An option file provides responses to all the questions that the platform poses. The format is described in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide), but it's basically a text file, 1 entry per line, with support for comments. It can be fed to the platform using the -o option.

The platform will also accept XML files containing a full description of a run. The schema for these files is given in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide). These files are also produced by the platform after every simulation, and placed in xml/<simulation ID>.xml. So if you want to rerun a simulation, you can point the platform to the XML file from that simulation using the -x option. Note that a new simulation ID will be assigned to the run, so there is no risk of overwriting previous simulation results.

### Available Options

To get a list of the current available options, run run_bbp.py with the -h flag.

```
 $ ./run_bbp.py --help
 Usage: run_bbp.py [options]
 
 Options:
  -h, --help                                show this help message and exit
  -x XML_FILE, --xml-file=XML_FILE          Run using XML description of workflow
  -s SIM_ID, --sim-id=SIM_ID                Force a sim id
  -o OPTFILE, --option-file=OPTFILE         File containing responses to interactive platform prompts
  -v, --version                             Broadband platform version
  -g, --generate-only                       Generates the XML description but does not run the platform
  -l LOG_FILE, --log=LOG_FILE               Directs output to a file, use to run BBP in background
  -m, --no-xml                              Do not generate xml
  -r RESUME_MODULE, --resume=RESUME_MODULE  Resume workflow from a certain module
  -e END_MODULE, --end=END_MODULE           End workflow after a certain module
  --expert                                  Turns on expert mode
```

### Validation Simulations

To run a validation simulation, go to the data/run directory and run run_bbp.py. The platform will ask you a series of questions. Answer 'y' to "Do you want to perform a validation run?"

```
 $ cd $BBP_DATA_DIR/run
 $ run_bbp.py
 Welcome to the SCEC Broadband Platform version 17.3.0.
 ================================================================================
 
 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event
 
 Do you want to perform a validation simulation (y/n)? y
 ================================================================================
 
 Please select a validation event from the list below:
 
 (1) Northridge
 (2) Loma Prieta
 ?
 ...
```

No input files are required by the user. However, you may wish to customize the validation simulation by selecting an alternate source description (src file) or a reduced station list to speed up the computations. You can put your own source description andor station list into the run directory (the format is described in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide)) or you can tell the platform where each file is located by using an absolute path. Note that any stations which do not have observed seismograms will not be included in the automatically generated goodness-of-fit comparison.

In addition to the low-frequency modules which compute seismograms using 1D Green's Tensors, certain methods allow validation events to be also run using precomputed seismograms to supply the low-frequency.

### User-defined Simulations

To run a user-defined simulation, two input files are required, a source description (src file) and a station list (stl file). For certain methods, the source description can either be in SRF format or a simple source description (the format is described in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide)). Other methods will always require a simple source description (src file). To run a user-defined simulation, run run_bbp.py.

```
 $ cd $BBP_DATA_DIR/run
 $ run_bbp.py
 Welcome to the SCEC Broadband Platform version 17.3.0.
 ================================================================================
 
 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event
 
 Do you want to perform a validation simulation (y/n)? n
 ================================================================================
 
 The Broadband Platform provides the following velocity models, which also include
 several method-specific and region-specific parameters.
 
 Please select a velocity model (either number or name is ok):
 
 (1) CEUS1k
 (2) CentralJapan
 (3) Canada1k
 (4) CentralUS
 (5) NOCAL
 (6) LABasin
 (7) WesternJapan
 (8) Mojave
 ?
 ...
```

You may then choose the method you would like to run:

```
 The Broadband Platform includes several scientific methods that can be
 used to calculate synthetic seismograms.
 
 Choose a Method to use in this Broadband scenario simulation:
 (1) GP (Graves & Pitarka)
 (2) UCSB
 (3) SDSU
 (4) EXSIM
 (5) CSM (Composite Source Model) - Beta Version
 (6) Song
 (7) Irikura Recipe Method 1
 ?
```

### Logging

During the run, log files will be produced in `logs/<simulation ID>/<simulation ID>.<module name>.log`.  If the platform fails, this is a good place to look to determine the error. Additionally, any fatal errors will be recorded in fatal_error.log.

Metadata capturing all the executable calls is located in `tmpdata/<simulation ID>/metadata.txt` for careful tracing of exactly what was called. Both the log files and metadata can be useful if troubleshooting an issue.

### Data Products

The platform produces a variety of data products. All data products are located in `outdata/<simulation ID>`. The ouput folder has a `index-<simulation ID>.html` file with a listing of all the data products in the output folder. You can open the index file on your browser window and then click through all data products. Image files can be displayed by clicking on the link to the file in the index.html page. The .bbp and .list files will be displayed as text files in the browser. On Mac OS X, you can see these data products by opening the outdata folder in Finder and double clicking on the specific file. On most Linux systems, you can show images using display:

```
 $ display <PNG file>
```

Make sure you have X11 forwarding enabled.

#### Station Map

To help visualize the stations in relationship to the fault, the platform produces a PNG file displaying station locations with red circles and the fault plane with a black line. You can find this file in `outdata/<simulation ID>/station_map.png`.

The Broadband Platform also produces a station_map.kml file that can be opened in Google Earth and will show the fault location and stations with their labels. This file is useful to identify the individual stations on the map.

#### Seismograms

When running the Broadband Platform, you have the option to generate plots of velocity and acceleration seismograms for each station. Plots of these files can be found in `outdata/<simulation ID>/<simulation ID>.<station>_<velocity or acceleration>_seis.png`.

The raw seismogram data is available in `outdata/<simulation ID>/<simulation ID>.<station>.vel.bbp` (velocity) and `outdata/<simulation ID>.<station>.acc.bbp` (acceleration).  Its format is described in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide).

#### Response Spectra

The RotD50 code, run at the end of each simulation, calculates the response spectra for each station. The raw RotD50 data is located at

```
 outdata/<simulation ID>/<simulation ID>.<station>.rd50
```

in the format described in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide).

#### Goodness-of-fit

If you run a goodness-of-fit module, several additional data products are produced. Two goodness-of-fit modules are available in Broadband.

* GP Goodness-of-fit: The goodness-of-fit comparison is performed by comparing the response spectra of a set of calculated seismograms to seismograms from another simulation or observed seismograms. For each station involved in the comparison, a plot comparing the response spectra calculated by the RotD50 module can be found in `outdata/<simulation ID>/<comparison label>_<simulation ID>/<station>rotd50.png`. A plot showing the seismograms on top and bottom can be found at `outdata/<comparison label>_<simulation ID>_<station>_overlay.png`. The seismogram comparison plot also includes a comparison of the arias duration. Finally, a goodness-of-fit plot with data from RotD50 is generated. It can can be found at `gof-<comparison label>-<simulation ID>_r0-<cutoff distance>-rotd50.png`.

   Note that at least 3 stations must be run for goodness-of-fit to be valid. If fewer than 3 stations are run, no goodness of fit calculation will be performed.

* SDSU Goodness-of-fit: The SDSU goodness-of-fit module can compute several GoF measures. The goodness-of-fit comparison is performed by comparing the synthetic seismogram to the observed seismogram for each of the sites. For each GoF measure computed, a summary file labeled `GOF_<measure>.list` is generated with the GoF values for each of the sites. A weighted sum of GoF measures called site-fit is available in gof_summary.txt file. If the simulation was run with three or more sites in the station list, GoF map plots will be available in the output directory as `<simulation ID>_<Gof measure>_map.png`. If PGA or PGV GoF values are computed a set of overlay plots showing the observed and synthetic seismograms on top and bottom will be available as `outdata/gof_plots/<simulation ID>_<station>_match_<format>_overlay.png` files.

#### Rupture files

When a user-defined event is simulated, the user has the option to run a rupture generator. This generator produces an SRF file, found in `outdata/<simulation ID>/*.srf`.  This file can be put in the run directory and used in future runs.  Additionally, the platform produces a plot of the cumulative slip on the fault surface, `outdata/<simulation ID>/<SRF prefix>.png`.

### Examples

Below are some examples that you can try using the sample files in the examples directory. You will need to have the Northridge Validation Package installed to use this example, which also requires that the LABasin Green's Functions package be installed in the Broadband Platform. Make sure all the tests pass before you try this. 

You should be in the run directory when you start these examples. If you set the aliases defined above, type:

```
$ run
```

Otherwise, type a path to the run directory. The Run directory provides a collection point, a staging area, when you are assembling the input files to be used in a simulation:

```
 $ cd $BBP_DATA_DIR/run
```

#### Running a Northridge Full Validation Simulation

You don't need to move any files for this. Notice that we will be using the EXSIM method, which will generate results quickly.

```
 $ run_bbp.py
 Welcome to the SCEC Broadband Platform version 17.3.0
 ================================================================================
 
 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event
 
 Do you want to perform a validation simulation (y/n)? y
 ================================================================================
 
 Please select a validation event from the list below:
 
 (1) Loma Prieta
 (2) Northridge
 ? 2
 ================================================================================
 
 The Broadband Platform includes several scientific methods that can be
 used to calculate synthetic seismograms.
 
 Choose a Method to use in this Broadband validation simulation:
 (1) GP (Graves & Pitarka)
 (2) UCSB
 (3) SDSU
 (4) EXSIM
 (5) CSM (Composite Source Model) - Beta Version
 (6) Song
 (7) Irikura Recipe Method 1
 ? 4
 ================================================================================
 SRC file: /home/sarah/bbp/bbp_val/Northridge/common/nr_v14_02_1.src
 ================================================================================
 STL file: /home/sarah/bbp/bbp_val/Northridge/common/nr_stats-ver_9_2013.stl
 Running ExSim
 ...
```

This simulation should take about 30 minutes (as opposed to about 12+ hours if using some of the other methods). Once it's complete, the platform will tell you:

```
 You can find results in $BBP_DATA_DIR/outdata/<simulation ID>
```

In that directory you will find:

* HTML Directory Listing (`index-<simulation ID>.html`)
* Velocity seismograms (`<simulation ID>.<station>.vel.bbp`)
* Acceleration seismograms (`<simulation ID>.<station>.acc.bbp`)
* Plots of velocity seismograms (`<simulation ID>.<station>_velocity_seis.png`)
* Plots of acceleration seismograms (`<simulation ID>.<station>_acceleration_seis.png`)
* Response spectra RotD50 files (`<simulation ID>.<station>.rd50`)
* Plots comparing simulated and observed seismograms (`Northridge_<simulation ID>_<station>_overlay.png`)
* Plots comparing simulated and observed RotD50 response spectra (`Northridge_<simulation ID>_<station>_rotd50.png`)
* Overall RotD50 goodness-of-fit plots (`gof-Northridge-<simulation ID>_r0-120-rd50.png`)

#### Sample Northridge Validation simulation, custom stations

Validation runs can take a long time. The time needed to generate each low-frequency seismogram will generaly increase with the magnitude of the event, taking more than 1 hour per station for earthquakes with Mw > 7.2. Sometimes you might want to run with a reduced station list so the simulation will run faster.

Copy the station list from `$BBP_DIR/doc/examples/northridge_3_stations` into the run directory, which should be `$BBP_DATA_DIR/run`. Take a look at the format of the station file:

```
 $ more northridge_3_stations.stl
 #Required:  lon, lat, station name, Vs30
 #Optional: low freq corner, high freq corner
 #lon		lat	station	Vs30	LF corn	HF corn
 -118.6417       34.5640 cast    450     0.120   111.11
 -118.4180       34.0628 lacn    278     0.140   111.11
 -118.8811       34.2886 moor    405     0.160   111.11
```

Now, run the platform, using a different station list. To do this, you will need to use the Broadband Platform "expert" mode. See below:

```
 $ run_bbp.py --expert
 Welcome to the SCEC Broadband Platform version 17.3.0.
 ================================================================================
 
 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event
 
 Do you want to perform a validation simulation (y/n)? y
 ================================================================================
 
 Please select a validation event from the list below:
 
 (1) Loma Prieta
 (2) Northridge
 ? 2
 ================================================================================
 
 The Broadband Platform includes several scientific methods that can be used to
 calculate synthetic seismograms.
 
 Choose a Method to use in this Broadband validation simulation:
 (1) GP (Graves & Pitarka)
 (2) GP Seis (using precomp seismograms)
 (3) UCSB
 (4) SDSU
 (5) SDSU Seis (using precomp seismograms)
 (6) EXSIM
 (7) CSM (Composite Source Model) - Beta Version
 (8) Song
 (9) Irikura Recipe Method 1
 ? 1
 ================================================================================
 
 When starting a simulation from a source description (SRC) file, the Broadband Platform
 workflow should include a rupture generator. Answer 'yes' here unless providing
 a complex Standard Rupture Format (SRF) file.
 
 Do you want to run the rupture generator (y/n)? y
 ================================================================================
  
 Each validation package includes a default source description (SRC) file for a historical event.
 Would you like to provide a different file instead of the default file provided?
 Answer 'no' here if you would like to use the standard source file for this event.
 
 Do you want to provide a custom source file (y/n)? n
 ================================================================================
 SRC file: /home/sarah/bbp/bbp_val/Northridge/common/nr_v14_02_1.src
 ================================================================================
 
 Station Selection
 =================
 Would you like to:
    (1) generate seismograms for all stations in the validation package
        OR
    (2) provide a custom list with a subset of the stations
 ? 2
 Do you want to
    (1) select a BBP station list in /home/sarah/bbb/bbp_data/run
        OR
    (2) enter the path of a BBP station list file
 ? 1
```

You will see a list of station list files that you have in your run directory. For this example, we only have the `northridge_3_stations.stl` file that we copied from the examples directory. You could have more files in the run directory if you have already tried other unit and acceptance tests. For this example, just select the `northridge_3_stations.stl` file.

```
 Here are the BBP station list files in the run directory.
 Please select one: 
 
 (1) northridge_3_stations.stl
 ? 1
 ================================================================================
 STL file: /home/sarah/bbp/bbp_data/run/northridge_3_sta.stl
 ================================================================================
 
 Site Response
 =============
 Running a site response module is an optional step while running a
 Broadband Platform simulation. It requires a station list file containing
 the Vs30 values for each station location.
 
 Do you want to run the site response module (y/n)? n
 ================================================================================
 Do you want to generate velocity seismograms' plots (y/n)? y
 ================================================================================
 Do you want to generate acceleration seismograms' plots (y/n)? y
 ================================================================================
 
 The Broadband Platform can generate comparison plots of the validation data against
 GMPEs to show how GMPEs match the recorded data for a certain event.
 
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
 Broadband Platform simulation. It creates a comparison plot showing how
 well the calculated seismograms fit recorded data.
 
 Do you want to run a goodness-of-fit module (y/n)? y
 ================================================================================
 
 Users can optionally select a Goodness of Fit module to plot a comparison of how well
 the simulated seismograms match the recorded data in a historical event.
 
 Choose a Goodness of Fit (GOF) Module:
 (1) GP
 (2) SDSU
 (3) Both
 ? 1
 ================================================================================
 
 Additional Metrics
 ==================
 Calculating additional metrics is an optional step on the Broadband Platform. It creates
 additional plots and data files that can be used to study the simulation.
 
 Do you want to calculate additional metrics (y/n)? n
```

Again, when the run completes, in about 15 minutes, you can find results in the output directory. You'll notice far fewer files, as only 3 stations were run instead of 133. The goodness-of-fit plots won't look very good - more stations are really needed to get an accurate plot. In the `$BBP_DIR/doc/examples/northridge_3_stations/results` directory you will find pre-calculated results for a few of the methods (GP and UCSB). You should obtain results equivalent to those.

#### Sample User-defined Southern California simulation with source description ==

Next let's try running a user-defined event. Copy the files `$BBP_DIR/examples/scenario_1_station/nr_one_stat.stl` and `$BBP_DIR/examples/scenario_1_station/user_eq.src` to the `$BBP_DATA_DIR/run` directory. user_eq.src is a simple source description. Its format is outlined in [File Formats](https://github.com/SCECcode/BBP/wiki/File-Format-Guide).

```
 $ run_bbp.py 
 Welcome to the SCEC Broadband Platform version 17.3.0.
 ================================================================================
 
 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event
 
 Do you want to perform a validation simulation (y/n)? n
 ================================================================================
 
 The Broadband Platform provides the following velocity models, which also include several
 method-specific and region-specific parameters.
 
 Please select a velocity model (either number or name is ok):
 
 (1) Mojave
 (2) LABasin
 (3) CentralJapan
 (4) LOMAP
 (5) WesternJapan
 ? 2
```

Since we are running a Southern California event, we select the LABasin velocity model, as it is the one that most closely matches the location of our rupture and stations.

```
 ================================================================================
 
 The Broadband Platform includes several scientific methods that can be used
 to calculate synthetic seismograms.
 
 Choose a Method to use in this Broadband scenario simulation:
 (1) GP (Graves & Pitarka)
 (2) UCSB
 (3) SDSU
 (4) EXSIM
 (5) CSM (Composite Source Model) - Beta Version
 (6) Song
 (7) Irikura Recipe Method 1
 ? 1
 
 ================================================================================
 
 The source description (SRC) file contains a description of the hypothetical
 (or scenario) earthquake, including information like location, geometry,
 magnitude, and mechanism.
 
 Do you want to
    (1) select a source description in /home/sarah/bbp/bbp_data/run
        OR
    (2) enter the path of a source description file
 ? 1
 Here are the source description files in the run directory.
 Please select one: 
 
 (1) user_eq.src
 ? 1
 Do you want to
    (1) select a BBP station list in /home/sarah/bbp/bbp_data/run
        OR
    (2) enter the path of a BBP station list file
 ? 1
 
 Here are the BBP station list files in the run directory.
 Please select one: 
 
 (1) nr_one_stat.stl
 (2) nr_three_stat.stl
 ? 1
 Running Genslip
 ...
```

Since this run only includes one station, it will run in about 5 minutes. In the output directory you'll notice there are no goodness-of-fit or files, since this is not a validation simulation. However, there is also a map file (station_map.png and station_map.kml), showing the fault plane and the stations, and a plot of the rupture slip (user_eq.png). The SRF generated by the rupture generator is in user_eq.srf; this file could be used in future runs. The filenames of the rupture slip plot and SRF are taken from the rupture description filename.