The Broadband Platform is open-source software that is made available under the terms of the Apache License, Version 2.0. A copy of the License is provided by the Apache Software Foundation (http://www.apache.org/licenses/LICENSE-2.0) and can also be found with the software in the 'LICENSE' file.

### Methods

To perform a simulation, a user selects a method, including optional modules. Currently, the Broadband Platform supports the following methods. Each method consists in a combination of modules (section below). A brief summary of each method is provided in the PDFs linked below, which include references pointing to more detailed publications:

* [Graves & Pitarka (GP)](pdfs/BBP-GP-2019-04-02.pdf)
* [SDSU](pdfs/BBP-SDSU_2019-04-02.pdf)
* [UCSB](pdfs/BBP-UCSB-2019-04-09.pdf)
* [EXSIM](pdfs/BBP-ExSIM-2019-04-09.pdf)
* [SONG](pdfs/BBP-Song-2019-04-03.pdf)
* [Irikura Recipe Method 1](pdfs/BBP-Irikura_Recipe_Method_1-2019-04-03.pdf)
* [Irikura Recipe Method 2](pdfs/BBP-Irikura_Recipe_Method_2-2019-04-11.pdf)
* Composite Source Model (CSM - Under development)

### Modules

The Broadband Platform consists of a series of modules. There are two main types of modules, science modules and utility modules. Science modules are those for which the platform has multiple implementations, provided by different coding research groups. The Science modules are combined to define a method (section above). Utility modules only have 1 implementation and are used by all simulations. A schematic of the available modules and their flow relationships is shown below in the following sections.

#### Science Modules

All simulations must include a module that creates synthetic seismograms. In some methods, there are separate low-frequency and high-frequency modules, while in other methods, these two steps are done by a single module. Rupture generation and site response are optional science modules. Users may select the following different implementations of each of these modules:

* Rupture generation: GP, UCSB, SONG, Irikura Recipe Method 1
* Low-frequency: GP, UCSB, EXSIM
* High-frequency: GP, UCSB, SDSU, EXSIM, Irikura Recipe Method 2
* Site response: GP

#### Post-processing Modules

A spectral response post-processing module is automatically run after the seismogram synthesis and the optional site response module are completed. For validation simulations, where simulated seismograms are compared against recorded data, users may select an optional goodness-of-fit (GoF) utility module to run at the conclusion of the simulation. The Broadband Platform currently supports both the GP and SDSU GoF modules and the users can select to run one, or both (or none) of these modules. Additionally, in validation simulations, users may also select to calculate the following optional metrics:

* RZZ2015 Metrics
* Fourier Amplitude Spectrum (FAS)
* Afshari and Stewart 2016 GMPE
* RotD100
* Anderson GoF 2004
