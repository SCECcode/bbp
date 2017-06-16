### Methods

To perform a simulation, a user selects a method, including optional modules. Currently, the Broadband Platform supports the following methods:

* Graves & Pitarka (GP)
* SDSU
* UCSB
* Composite Source Model (CSM - Under development)
* EXSIM
* SONG
* Irikura Recipe Method 1

Our implementation is open-source and made available to users for the purposes of academic and technical research. We are not responsible for the usage of our codes for any other application.

### Modules

The Broadband Platform consists of a series of modules. There are two main types of modules, science modules and utility modules. Science modules are those for which the platform has multiple implementations, provided by different coding research groups. Utility modules only have 1 implementation and are used by all simulations. A schematic of the available modules and their flow relationships is shown below in the following sections.

#### Science Modules

All simulations must include a module that creates synthetic seismograms. In some methods, there are separate low-frequency and high-frequency modules, while in other methods, these two steps are done by a single module. Rupture generation and site response are optional science modules. Users may select the following different implementations of each of these modules:

* Rupture generation: GP, UCSB, SONG, Irikura Recipe Method 1
* Low-frequency: GP, UCSB, CSM, EXSIM
* High-frequency: GP, UCSB, SDSU, CSM, EXSIM
* Site response: GP, UCSB

#### Post-processing Modules

A spectral response post-processing module is automatically run after the seismogram synthesis and the optional site response module are completed. For validation simulations, where simulated seismograms are compared against recorded data, users may select an optional goodness-of-fit (GoF) utility module to run at the conclusion of the simulation. The Broadband Platform currently supports both the GP and SDSU GoF modules and the users can select to run one, or both (or none) of these modules. Additionally, in validation simulations, users may also select to calculate the following optional metrics:

* RZZ2015 Metrics
* Fourier Amplitude Spectrum (FAS)
* Afshari and Stewart 2016 GMPE
* RotD100
* Anderson GoF 2004