# Step2 - Computing the overlap matrices and nonadiabatic couplings in Kohn-Sham and TD-DFT level of theory

In order to compute the nonadiabatic couplings (NACs) in Kohn-Sham (KS) and TD-DFT level of theory we need follow these steps:

## 1. CP2K input for electronic structure calculations

`cp2k_input_template.inp`

The input template to compute the electronic structure with the `RUN_TYPE ENERGY`. The following input should be able to compute the TD-DFT calculations, producing the cube files
and producing the PDOS files (if specified by user). We recommend the use of the above input template and if the user needs to change the input based on his/her needs it is better to perform the changes on this input file. One should make sure for production of the cube files, performing the TD-DFT calculations. 

Other required files for running the CP2K input file are basis set and pseudopotential files or any other files required to run the calculations. One should specify the full path to these files in the `cp2k_input_template.inp` file.

For `TDDFPT` section the number of excited states must be specified. The higher number of excited states needs higher computational cost, and also one needs to specify the cube files to be printed for higher number of KS states. In this section the `RESTART` is required to be set to `.TRUE.`. Also, the `WFN_RESTART_FILE_NAME` in this section should exist.

For `BASIS_SET_FILE_NAME` and `POTENTIAL_FILE_NAME` keyword the user should specify the full path to these files. Also, `WFN_RESTRAT_FILE_NAME` should exist in the input. This will increase the speed of calculations.

For `&MO_CUBES` section the number of occupied and unoccupied orbitals must be specified. 

## 2. CP2K input for electronic structure calculations
