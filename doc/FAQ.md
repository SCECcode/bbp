# Frequently Asked Questions (FAQ) for the SCEC Broadband Platform

* What it the most recent version of the Broadband Platform?
  * The most recent version of the Broadband Platform software is posted GitHub's releases page. We expect to make new releases of the Broadband Platform every 3-6 months.

* Which version of the Broadband Platform should I use?
  * Several previous versions of the broadband platform have been released by SCEC. New versions of the platform incorporate software updates, improvements, and bug fixes. For this reason, we recommend all new users should use the most recent version of the platform, as posted on the Broadband Platform wiki page. We recommend users migrate their work to the most recent version of the broadband platform at the first opportunity. Previous versions of the platform are still available primarily to support reproducibility of earlier results.

* How can I simulate earthquakes in California, Japan and Eastern North America?
  * The newest Broadband Platform distribution supports simulations for several 1D velocity models. These 1D velocity models were defined to represent structure in parts of California including Los Angeles Basin, Mojave Desert, Central California Coast, Eastern North America, including Virginia and Eastern Canada, and Western and Central Japan. Broadband Platform users can select a region, and simulate earthquakes in that region for one or more of the ground motion modeling software installed in the platform.

* What region should I use for an earthquake in California?
   * The Broadband Platform currently has 3 simulation regions for California: LA Basin, Mojave, and Northern California. Please take a look at this [map](http://hypocenter.usc.edu/research/bbp/california_gfs_16_5.pdf) for a suggestion of what region to use for simulations in California. Please note that these are only approximate and somewhat subjective, but basically, for ruptures located in the Southern California, users should use the LA Basin region. The Mojave region should be used for faults in most of the inland part of Southern California, and the Northern California region should be used for ruptures located in the remainder of the State.

* How can I simulate earthquakes outside California, Eastern North America, and Japan?
  * If users want to simulate an earthquake outside of California, Eastern North America, or Japan, they must select a velocity profile for their region of interest. In future releases of the broadband platform, it may be possible for the user to define a 1D velocity profile, generate Greens Functions for that profile, and then simulate earthquakes using those Greens Functions. Currently, the only platform does not include software to generate new Greens functions, for most simulation methods (EXSIM is an exception to this rule).  So, to use most methods, users should examine the 1D velocity profiles for each of the supported 1D models currently supported, and select the velocity model most like their region of interest. Then, then can run simulations using existing Greens functions. Given a 1D velocity model, simulations can be run with faults can be embedded within this 1D region, and stations can be distributed around the fault, and the ground motion simulations results can be representative of other regions.

* I try to run the tests of UniTests.py. The result tells me that ImportError: No module named scipy. I wonder what wrong it is.
  * The error you are seeing relates to the python software that you are using. The Broadband platform using standard python libraries. It also requires several "optional" python libraries that must be added to most python installations. The error you are getting indicates that you need to install one (or more) additional python libraries. These additional libraries are free, but how you install them will depend on the system on which you are installing the broadband platform.

    We have three suggestions on how to solve this problem.

    First, possibly the easiest way to obtain a working version of the broadband platform is to use the Virtual box installation. In this approach, you install a free software tool, called virtual box, on your computer. Then you retrieve a Broadband platform virtual box image from SCEC that contains all necessary Broadband platform software including all python libraries. You load the broadband platform image into the virtual box, and you can run the broadband platform without adding any additional software. More details on using this approach are provided here:

    * http://scec.usc.edu/scecpedia/BBP_16.5.0_Virtual_Box_Image

    Second, you can install the necessary additional python libraries onto your current computer. Then the broadband platform should run. How to install the needed libraries, such as the scipy library mentioned in your error message will depend on the python software you are using. We recommend using the anaconda python library because it includes most of the required python libraries. We describe how to install anaconda libraries (for mac) in this entry. Installing an using anaconda python should be similar for other systems also. Here is a description from our website:

    * http://scec.usc.edu/scecpedia/BBP_on_OS_X_Guide_16_5

    The BBP v16.5.0 requires Python v2.7+ and several optional python libraries. It is possible to collect all necessary Python libraries, but we recommend using a Python distribution that includes nearly all the required libraries. We recommend use of the Continuum Anaconda Python distribution (free). Go to the Continuum Anaconda download site, and download the Python 2.7 command line distribution tool (free). Continuum Anaconda Python for Mac Once this is downloaded, run the installer by typing:

    ```% bash Anaconda2-4.0.0-MacOSX-x86_64.sh```

    It will ask you to agree to the license, and then install a current version of Python 2.7 and most of the Python libraries required by the BBP v16.5.0. Anaconda Python can be installed in your user directory and does not require administrator privileges. The end of the installation process will update your .bash_profile to put the Anaconda Python into your default path. We tested using the Anaconda version given above, and it worked with BBP v16.5.0 

    Update the PyProj library in Anaconda Python

    One Python package that must be updated to in the pyproj library. Use the Anaconda update tool. Open a terminal window and type:

     ```% conda install pyproj```

    The conda installer program will run, and it will identify what needs to be installed, or updated. When it asks for approval to proceed, enter yes. Pyproj will be added/updated to your Mac.

    Your third option is to use the python updating utilities currently installed on your computer, and add the scipy library. One python installer is called pip.

    ```$ pip install scipy```

    If you use this approach, you may find that you need to add other python modules (in addition to scipy) such as matplotlib, pyproj, and basemap. However, once you add all the required python modules, the Broadband Platform software should run. 

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

    You are receiving these warnings are because you are using a new version of one of the Python modules. These are warning and not errors, and they will not affect the results you get using the Broadband Platform. We will expect to change the Broadband software to eliminate these warnings in the next Broadband release.