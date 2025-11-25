### Modules

The Broadband Platform consists of a series of modules. There are two main types of modules, science modules and utility modules. Science modules are those for which the platform has multiple implementations, provided by different coding research groups. The Science modules are combined to define a method (see sub-section below). Utility modules only have 1 implementation and are used by all simulations. A schematic of the available modules and their flow relationships is shown below in the following sections.

#### Science Modules

All simulations must include some description of the source and a process by which to compute synthetic seismograms. We refer to the source description or modeling module as rupture generators (RG). In some methods, the seismogram computation is performed using separate low-frequency (LF) and high-frequency (HF) modules, while in other methods, these two steps are performed within a single module. Site response (SR) is an optional science module used to modify the seismograms to account for site conditions different that those computed based on reference velocity models. The science modules currently implemented in the BBP are:

* Rupture generation (RG): GP, UCSB, EXSIM, SONG, Irikura Recipe
* Low-frequency (LF): GP, UCSB, EXSIM
* High-frequency (HF): GP, UCSB, SDSU, EXSIM, Irikura Recipe Method 2
* Site response (SR): GP, PySeismoSoil

### Methods

The BBP was originally designed to provide users with a wide flexibility on possible module (section above) combinations. However, the BBP evolved into suites of preferred module combinations, which we refer to as "methods". Currently, BBP users can compute simulations from any of 7 distinct methods (listed below with their associated modules). A brief summary of each method is provided in the PDFs linked below, which includes references pointing to more detailed publications:

* [Graves & Pitarka (GP): uses RG, LF, and HF from GP](pdfs/GP_method-20220216.pdf)
* [SDSU: uses RG and LF from GP and HF from SDSU](pdfs/SDSU_release_2022.pdf)
* [UCSB: uses RG, LF, and HF from UCSB](pdfs/UCSB_21.3_c.pdf)
* [EXSIM: uses RG, LF, and HF form EXSIM](pdfs/EXSIM_20220309_V2.pdf)
* [SONG: used RG from SONG; LF and HF from GP](pdfs/BBB-Song-2022-03-02.pdf)
* [Irikura Recipe Method 1: uses the RG from the Irikura Recipe; LF and HF from GP](pdfs/Irikura_Recipe_Method_1.2022.pdf)
* [Irikura Recipe Method 2: uses the RG and HF from the Irikura Recipe; LF from GP](pdfs/Irikura2_2022v4.pdf)
* Composite Source Model (CSM - Under development)

There are currently two site response modules, GP and PySeismoSoil, which can be used by any of the simulation methods above, whenever site response is selected by the user.

#### Post-processing Modules

A spectral response post-processing module is automatically run after the seismogram synthesis and the optional site response module are completed. For validation simulations, where simulated seismograms are compared against recorded data, users may select an optional goodness-of-fit (GoF) utility module to run at the conclusion of the simulation. The Broadband Platform currently supports the GP, FAS, and SDSU GoF modules and the users can select to run one, two, all (or none) of these modules. Additionally, in validation simulations, users may also select to calculate the following optional metrics:

* RZZ2015 Metrics
* Afshari and Stewart 2016 GMPE
* RotD100
* Anderson GoF 2004

### Software License
The SCEC-developed portions of the Broadband platform software is distributed under the terms of the BSD 3-Clause open-source license. Please see the LICENSE.txt file for more information. Individual models codes may be offered under their own open-source software licenses, please look for a LICENSE.txt file under specific sub folders in the 'src' directory for module-specific license information.
