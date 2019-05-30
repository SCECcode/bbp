### Modules

The Broadband Platform consists of a series of modules. There are two main types of modules, science modules and utility modules. Science modules are those for which the platform has multiple implementations, provided by different coding research groups. The Science modules are combined to define a method (see sub-section below). Utility modules only have 1 implementation and are used by all simulations. A schematic of the available modules and their flow relationships is shown below in the following sections.

#### Science Modules

All simulations must include some description of the source and a process by which to compute synthetic seismograms. We refer to the source description or modeling module as rupture generators (RG). In some methods, the seismogram computation is performed using separate low-frequency (LF) and high-frequency (HF) modules, while in other methods, these two steps are performed within a single module. Site response (SR) is an optional science module used to modify the seismograms to account for site conditions different that those computed based on reference velocity models. The science modules currently implemented in the BBP are:

* Rupture generation (RG): GP, UCSB, EXSIM, SONG, Irikura Recipe
* Low-frequency (LF): GP, UCSB, EXSIM
* High-frequency (HF): GP, UCSB, SDSU, EXSIM, Irikura Recipe
* Site response (SR): GP

### Methods

The BBP was originally designed to provide users with a wide flexibility on possible module (section above) combinations. However, the BBP evolved into suites of preferred module combinations, which we refer to as "methods". Currently, BBP users can compute simulations from any of 7 distinct methods (listed below with their associated modules). A brief summary of each method is provided in the PDFs linked below, which includes references pointing to more detailed publications:

* [Graves & Pitarka (GP): uses RG, LF, and HF from GP](pdfs/BBP-GP-2019-04-02.pdf)
* [SDSU: uses RG and LF from GP and HF from SDSU](pdfs/BBP-SDSU_2019-04-02.pdf)
* [UCSB: uses RG, LF, and HF from UCSB](pdfs/BBP-UCSB-2019-04-09.pdf)
* [EXSIM: uses RG, LF, and HF form EXSIM](pdfs/BBP-ExSIM-2019-04-09.pdf)
* [SONG: used RG from SONG; LF and HF from GP](pdfs/BBP-Song-2019-04-03.pdf)
* [Irikura Recipe Method 1: uses the RG from the Irikura Recipe; LF and HF from GP](pdfs/BBP-Irikura_Recipe_Method_1-2019-04-03.pdf)
* [Irikura Recipe Method 2: uses the RG and HF from the Irikura Recipe; LF from GP](pdfs/BBP-Irikura_Recipe_Method_2-2019-04-11.pdf)
* Composite Source Model (CSM - Under development)

There is currently one site response module (GP) and it is used for all the methods above, whenever site response is selected by the user.

#### Post-processing Modules

A spectral response post-processing module is automatically run after the seismogram synthesis and the optional site response module are completed. For validation simulations, where simulated seismograms are compared against recorded data, users may select an optional goodness-of-fit (GoF) utility module to run at the conclusion of the simulation. The Broadband Platform currently supports both the GP and SDSU GoF modules and the users can select to run one, or both (or none) of these modules. Additionally, in validation simulations, users may also select to calculate the following optional metrics:

* RZZ2015 Metrics
* Fourier Amplitude Spectrum (FAS)
* Afshari and Stewart 2016 GMPE
* RotD100
* Anderson GoF 2004

### Software License
The Broadband Platform is open-source software that is made available under the terms of the Apache License, Version 2.0. A copy of the License is provided by the Apache Software Foundation (http://www.apache.org/licenses/LICENSE-2.0) and can also be found with the software in the 'LICENSE' file.
