# Step2 - Computing the overlap matrices and nonadiabatic couplings in Kohn-Sham and TD-DFT level of theory

In order to compute the nonadiabatic couplings (NACs) in Kohn-Sham (KS) and TD-DFT level of theory we need follow these steps:

## 1. CP2K input for electronic structure calculations

`cp2k_input_template.inp`

This input is a template to compute the electronic structure with the `RUN_TYPE ENERGY`. The input should be able to compute the TD-DFT calculations, producing the cube files
and PDOS files. We recommend the use of the above input template and if the user needs to change the input based on his/her needs it is better to perform the changes on this input file. One should make sure about the production of the cube files and performing the TD-DFT calculations.

Other required files for running the CP2K input file are basis set and pseudopotential files or any other files required to run the calculations, such as `dftd3.dat`. The full path to these files in the `cp2k_input_template.inp` file shoud be specified.

For `TDDFPT` section the number of excited states are specified. Higher number of excited states needs higher computational cost, and also one needs to specify the cube files to be printed for higher number of KS states. For more information about TD-DFPT calculations please refer to the following links:
[CP2K paper](https://aip.scitation.org/doi/pdf/10.1063/5.0007045)
[Difference between TD-DFT and TD-DFPT](https://groups.google.com/g/cp2k/c/xj8udnSyeEI)

The keyword `RESTART` increase the speed of calculations, both for SCF and TD-DFPT calculations. Therefore, the `RESTART` is required to be set to `.TRUE.` in `TDDFPT` section. Also, the `WFN_RESTART_FILE_NAME` in this section should exist with a random `tdwfn` file name. The same is also needed in the `FORCE_EVAL` section. `WFN_RESTRAT_FILE_NAME` should exist in the input with an random `wfn` file name. 

In the `&MO_CUBES` section the number of occupied and unoccupied orbitals must be specified. This is dependent on the TD-DFPT calculations and the user have to make a good guess to make sure that the cube files of all the states in the excitation analysis of TD-DFPT calculations exist.

## 2. 
