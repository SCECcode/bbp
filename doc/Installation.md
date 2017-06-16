Installing the Broadband Platform involves obtaining a copy of the code and building the required executables. You can either download the tar.gz platform from BBP's GitHub releases page or check the code out of GitHub directly. Users can also obtain a Virtual Box Linux image with everything they need to run the Broadband Platform.

### Software Dependencies

The Broadband Platform has the following software dependencies:

* Python v2.7.9+ with
  * Matplotlib 1.4.3+
  * NumPy 1.9.2+
  * SciPy 0.14.1+
  * Pyproj 1.9.2+
* GNU compilers (gcc, gfortran) 4.5.1+

Please make sure they are installed in your computer before you continue with the BBP installation process.

### Easy Installation

Most users will want to use the easy installation procedure described below to install the Broadband Platform to their systems. We recommend that users that want to use the Broadband Platform download the tar.gz file available on GitHub's release page. For advanced users, who would like to make modifications to the BBP and contribute these modifications back to us, we recommend they clone our repository so that it is easier to track their changes to the software.

To install the Broadband Platform on your computer using the easy installation script, please follow these steps:

* Create a directory on your computer where you want all Broadband packages to be installed:

```
  $ cd /home/sarah
  $ mkdir bbp
  $ cd /home/sarah/bbp
```

* Download the .tar.gz file from GitHub into the recently-created directory
* Uncompress the downloaded file

```
  $ tar -xzf bbp-17.3.0.tar.gz
```

* Delete the tar.gz file as it will not be needed anymore (optional)

* Run the easy installation script located in the setup sub-directory. For members of the SCEC community, please use your SCEC account to [log into our wiki](https://scec.usc.edu/it/Broadband_Easy_Install_v17.3.0) to download another installation script to install additional validation packages for the BBP with a larger number of events.

```
  $ cd bbp-17.3.0/setup
  $ ./easy_install_bbp.sh
```

The easy installation script will create a number of other directories inside your top-level bbp directory (/home/sarah/bbp). It will compile the codes in the BBP source distribution, and it will then ask the user which BBP packages should be installed. Each package contains a region (e.g. LA Basin, Northern California, etc). The user should select 1 (Yes) to install the package, or 2 (No) to skip it. Please note that the LA Basin region is required for running the Unit and Acceptance tests provided with the Broadband Platform.

After asking the user about each of the 8 available regions, the easy install script will download and install of these regions. At the end of the installation process, the script will print a number of commands that should be inserted in the user's .bash_profile file. In the example above, you will see an output like:

```
=> All Done!

Please add the following lines to your bash_profile:

export BBP_DIR=/home/sarah/bbp/bbp-17.3.0/bbp
export BBP_GF_DIR=/home/sarah/bbp/bbp_gf
export BBP_VAL_DIR=/home/sarah/bbp/bbp_val
export PYTHONPATH=/home/sarah/bbp/bbp-17.3.0/bbp/comps
export BBP_DATA_DIR=/home/sarah/bbp/bbp_data
export PATH=/home/sarah/bbp/bbp-17.3.0/bbp/comps:/home/sarah/bbp/bbp-17.3.0/bbp/utils/batch:$PATH
ulimit -s unlimited
```

After installing the Broadband Platform on their systems, users should confirm the code is built correctly by [running Unit and Acceptance Tests](https://github.com/SCECcode/BBP/wiki/Running-Tests) before starting to use the code for research purposes.