# The SCEC Broadband Platform (BBP) Software

<a href="https://github.com/sceccode/bbp.git"><img src="https://github.com/sceccode/bbp/wiki/images/SRL_Cover_v8.png"></a>

[![Python](https://img.shields.io/badge/python-%3E%3D3.7-blue)](https://www.python.org)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/bbp)
[![bbp-ci Actions Status](https://github.com/SCECcode/bbp/workflows/bbp-ci/badge.svg)](https://github.com/SCECcode/bbp/actions)

## Description
The Southern California Earthquake Center (SCEC) Broadband Platform (BBP) is a software system that can generate 0-20+ Hz seismograms for historical and scenario earthquakes in California, Eastern North America, and Japan using several alternative computational methods.

## Table of Contents
1. [Software Documentation](https://github.com/SCECcode/bbp/wiki)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Support](#support)
5. [Citation](#citation)
6. [Contributing](#contributing)
7. [Credits](#credits)
8. [License](#license)

## Installation

BBP was developed to support ground simulations run on Linux servers and high-performance computing systems,
so it is designed to compile and run on Linux-based computers. Before installing BBP, they should be aware that it is
possible to run some versions of the BBP software without installing the software on a Linux computer.
Below we outline the two options for running the BBP software:

1. [BBP Docker Images](https://github.com/sceccode/bbp_docker) Users can run a version of BBP using Docker on their
local computers including laptops. The BBP Docker version contains an LA-Basin region velocity model, and three validation events in that regeion.
Users can install free Docker software on most computers (e.g. Linux, MacOS, Windows) then run a BBP Docker image in a terminal window on their computer.
2. [Installation Instructions for Linux Systems](https://github.com/SCECcode/bbp/wiki/Installation)
Advanced users that want to install many or all of the BBP simulation regions models, or that want to run large
parallel queries of the CVM models, should install the BBP software on a Linux system. BBP software is developed
on USC Center for Advanced Research Computing (CARC) Linux cluster.

## Usage

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

To run a validation simulation, go to the data/run directory and run run_bbp.py. The platform will ask you a series of questions. Answer 'y' to "Do you want to perform a validation run?" (the exact list of validation events users see when running the Broadband Platform depends on what validation events are installed on their computers):

```
 $ run_bbp.py
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
 (3) LOMAP
 (4) NR
 (5) Whittier
 ?
 ...
```

No input files are required by the user. However, you may wish to customize the validation simulation by selecting an alternate source description (src file) or a reduced station list to speed up the computations. You can put your own source description and/or station list into the run directory (the format is described in [File Formats](./File-Format-Guide)) or you can tell the platform where each file is located by using an absolute path. Note that any stations which do not have observed seismograms will not be included in the automatically generated goodness-of-fit comparison. To supply alternative source description and/or station list files, please run the Broadband Platform in 'expert' mode using the '--expert' command-line flag.

### User-defined Simulations

To run a user-defined simulation, two input files are required, a source description (src file) and a station list (stl file). A simple source description (src file) is always required, but, for certain methods, a source description in SRF format (the format is described in [File Formats](./File-Format-Guide) can be supplied as well and will be used for the seismogram computation modules in the Broadband Platform. To run a user-defined simulation, run run_bbp.py:

```
 $ run_bbp.py
 Welcome to the SCEC Broadband Platform version 22.4.0.
 ================================================================================

 Please select the Broadband Platform mode of operation:
    * Validation - Simulates a historical event
    * Scenario   - Runs a user-defined hypothetical event

 Do you want to perform a validation simulation (y/n)? n
 ================================================================================

 The Broadband Platform provides the following velocity models, which also include
 several method-specific and region-specific parameters.

 Please select a velocity model (either number or name is ok):

 (1) CentralJapan500
 (2) CentralCal500
 (3) NOCAL500
 (4) LABasin500
 (5) WesternJapan500
 (6) Mojave500
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
 (5) Song
 (6) Irikura Recipe Method 1 (Irikura1)
 (7) Irikura Recipe Method 2 (Irikura2)
 ?
```

## Support
Support for BBP is provided by that Southern California Earthquake Center (SCEC) Research Computing Group. This group supports several research software distributions including BBP. Users can report issues and feature requests using BBP's github-based issue tracking link below. Developers will also respond to emails sent to the SCEC software contact listed below.
1. [BBP Github Issue Tracker](https://github.com/SCECcode/bbp/issues)
2. Email Contact: software@scec.org

## Citation
References, citations, and acknowledgements help us obtain continued support for the development of the BBP software. If you use the BBP software in your research, please include the citation of the BBP paper in the references/bibliography section of your publication. This is more effective than you providing in-text acknowledgements.

* Preferred Reference: Maechling, P. J., F. Silva, S. Callaghan, and T. H. Jordan (2015). SCEC Broadband Platform: System Architecture and Software Implementation, Seismol. Res. Lett., 86, no. 1, doi: 10.1785/0220140125.

* Example Acknowlegement: We would like to acknowledge the use of the SCEC Broadband Platform Software (Maechling 2015) in this research.

Along with citing the BBP software, researchers should also cite the appropriate publication for any of the ground motion models they use in their research. Citations for individual ground motion methods are included in the [Credits.md](CREDITS.md) file in this repository.

## Contributing
We welcome contributions to the BBP software framework.
Geoscientists can register their ground motion models into BBP and software developers can
improve and extend the BBP software. An overview of the process for contributing seismic models or
software updates to the BBP Project is provided in the BBP [contribution guidelines](CONTRIBUTING.md).
BBP contributors agree to abide by the code of conduct found in our [Code of Conduct](CODE_OF_CONDUCT.md) guidelines.

## Credits
Development of BBP is a group effort. A list of developers that have contributed to the BBP Software framework
are listed in the [Credits.md](CREDITS.md) file in this repository.

## License
The SCEC-developed portions of the Broadband platform software is distributed under the BSD 3-Clause open-source license.
Please see the [LICENSE.txt](LICENSE.txt) file for more information. Individual models codes may be offered under their own open-source software licenses.
