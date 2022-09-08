# Frequently Asked Questions (FAQ) for the SCEC Broadband Platform

* What it the most recent version of the Broadband Platform?
  * The most recent version of the Broadband Platform software is posted GitHub's releases page. We expect to make new releases of the Broadband Platform every 6-12 months.

* Which version of the Broadband Platform should I use?
  * Several previous versions of the Broadband Platform have been released by SCEC. New versions of the Platform incorporate software updates, improvements, and bug fixes. For this reason, we recommend all new users should use the most recent version of the Platform, as posted on the Broadband Platform wiki page. We recommend users migrate their work to the most recent version of the Broadband Platform at the first opportunity. Previous versions of the Platform are still available primarily to support reproducibility of earlier results.

* Is the Broadband Platform compatible with Python 3?
  * The Broadband Platform was ported to Python 3 as of version 19.8. All Broadband Platform versions after 19.8 require Python 3 to run.

* How can I simulate earthquakes in California, Japan and Eastern North America?
  * The newest Broadband Platform distribution supports simulations for several 1D velocity models. These 1D velocity models were defined to represent structure in parts of California including Los Angeles Basin, Mojave Desert, Central California Coast, Southern Walker Lane (Southern Sierra Nevada region), Eastern North America, including Virginia and Eastern Canada, and Western and Central Japan. Broadband Platform users can select a region, and simulate earthquakes in that region for one or more of the ground motion modeling software installed in the platform.

* What region should I use for an earthquake in California?
   * The Broadband Platform currently has 5 simulation regions for California: LA Basin, Mojave,
   Central California, Southern Sierra Nevada, and Northern California. Please take a look at this [map](pdfs/california_gfs_22_4.pdf) for a suggestion of what region to use for simulations in California. Please note that these are only approximate and somewhat subjective, but basically, for ruptures located in the Southern California, users should use the LA Basin region. The Mojave region should be used for faults in most of the inland part of Southern California, the Central California region should be used for the central coast region of California, the Southern Sierra Nevada region should be used for inland portion of Central and Northern California, and the Northern California region should be used for ruptures located in the remainder of the State.

* How can I simulate earthquakes outside California, Eastern North America, and Japan?
  * Simulating earthquakes outside one of the currently supported regions requires creating a new 1D velocity profile and setting up a new simulation region. In addition to the 1D velocity profile for the region, the Broadband Platform requires the calculation of Greens' Functions and the configuration of several region-specific parameters so it is not as simple as just creating a new 1D velocity model and supplying that file to the BBP. At this time, setting up a new velocity model/Green's Functions/region parameters requires help from several of the science teams. Also the Broadband platform does not currently include software to generate new Greens functions needed for most simulation methods (EXSIM is an exception to this rule). In future releases of the Broadband Platform, it may be possible for the user to define a 1D velocity profile, generate Greens Functions for that profile, and then simulate earthquakes using those Greens Functions.

* Every time I run the Broadband Platform I always get the exact same results. What do I have to do to generate different results when I run the Broadband Platform?
  * The Broadband Platform was designed to generate reproducible results. Therefore, given the same inputs, we expect it to produce the exact outputs again and again. In order to generate different simulation variations (or realizations), users can run the Broadband Platform multiple times, but providing a different SEED value each time the BBP is run. The SEED value is included in the simple source description (SRC) file. So, to run 10 different simulation variations, users will need to create 10 copies of the SRC file, each with a different value for the SEED parameter. This seed will be used to initialize the random number generators in several of the BBP modules, resulting in different parameters to be used. For example, in the GP rupture generator, different SEED values will result in different slip distributions. Similarly, other modules use the SEED to randomize their own parameters.

* I try to run the tests of UnitTests.py. The result tells me that ImportError: No module named scipy. I wonder what wrong it is.
  * The error you are seeing relates to the python software that you are using. The Broadband Platform uses standard Python libraries. It also requires several "optional" Python libraries that must be added to most Python installations. The error you are getting indicates that you need to install one (or more) additional Python libraries. These additional libraries are free, but how you install them will depend on the system on which you are installing the broadband platform.

    Users can install the necessary additional Python libraries onto their computer. Then the Broadband Platform should run. How to install the needed libraries, such as the SciPy library mentioned in your error message will depend on the Python software you are using.

    One option that we recommend is using the open source [Anaconda](https://www.anaconda.com) Python distribution. It includes most of the required python libraries in one simple installation.

    The BBP v22.4.0 requires Python v3.6+ and several optional Python libraries. It is possible to collect all necessary Python libraries, but we recommend using a Python distribution that includes nearly all the required libraries. We recommend use of the [Anaconda Python distribution (free)](https://www.anaconda.com). Go to the Anaconda downloads page, and download the Python 3.6+ command line distribution tool (free). Once downloaded, the Anaconda Python distribution can be installed by running the downloaded script file.

    It will ask you to agree to the license, and then install a current version of Python 3.6+ and most of the Python libraries required by the BBP v22.4.0. Anaconda Python can be installed in your user directory and does not require administrator privileges. At the end of the installation process, the installation tool should update your shell initialization files to put the Anaconda Python installation into your default path.

    Update the PyProj library in Anaconda Python

    One Python package that must be updated to in the pyproj library. Use the Anaconda update tool. Open a terminal window and type:

    ```
    $ conda install pyproj
    ```

    The conda installer program will run, and it will identify what needs to be installed, or updated. When it asks for approval to proceed, enter yes. Pyproj will then be added/updated to your system.

    Another option that users have is to use the Python updating utilities currently installed on your computer, and add the scipy library. One Python installer is called pip.

    ```
    $ pip install scipy
    ```

    If you use this approach, you may find that you need to add other Python modules (in addition to scipy) such as matplotlib, pyproj, numba, and basemap. However, once you add all the required Python modules, the Broadband Platform software should run.

* When I run the UnitTests.py, I get a few warnings. I wonder if I can ignore these warnings, what do I do next?
  * When running BBP's unit tests, you will sometime see warning messages like these:
    ```
    test_rmg (test_rmg.TestRMG) ... /home/liyq/anaconda2/lib/python2.7/site-packages/numpy/core/fromnumeric.py:225: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
    return reshape(newshape, order=order)
    /home/liyq/anaconda2/lib/python2.7/site-packages/numpy/core/numeric.py:190: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
    a = empty(shape, dtype, order)
    /home/liyq/anaconda2/lib/python2.7/site-packages/numpy/lib/shape_base.py:873: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
    return c.reshape(shape_out)
    ```

    You are receiving these warnings are because you are using a newer version of one of the Python modules. These are warning and not errors, and they should not affect the results you get using the Broadband Platform. We will expect to change the Broadband software to eliminate these warnings in the next Broadband release.
