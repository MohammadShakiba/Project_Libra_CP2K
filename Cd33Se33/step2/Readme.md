# Step2 - Computing the overlap matrices and nonadiabatic couplings in Kohn-Sham and TD-DFT level of theory

In order to compute the nonadiabatic couplings (NACs) in Kohn-Sham (KS) and TD-DFT level of theory we need to have the following files:

## CP2K input for computing TD-DFT and

`cp2k_input_template.inp

The input template to compute the electronic structure with the `RUN_TYPE ENERGY`. The following input should be able to compute the TD-DFT calculations, producing the cube files
and poducing the PDOS files (if specified by user).
