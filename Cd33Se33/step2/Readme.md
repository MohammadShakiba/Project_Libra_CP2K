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

In the `&MO_CUBES` section the number of occupied and unoccupied orbitals must be specified. This is dependent on the TD-DFPT calculations and the user have to make a good guess to make sure that the cube files of all the states in the excitation analysis of TD-DFPT calculations exist. This guess can be obtained from running the calculations for 5-10 steps.

## 2. Bash file for running the calculations for one job

The standard sample bash file for submitting the calculations and running the Python code through `slurm` and `sbatch` is the `submit_template.slm`. First one should load all needed modules including the modules required for loading CP2K. This file contains the input variables required for calculations of the overlap matrices and NACs. We list the variables as follows:

`nprocs`: The number of processors used for calculations. Note that the same number of processors should be specified in the `#SBATCH --ntasks-per-node` above.

`cp2k_exe`: The executable CP2K or the full path to executable CP2K folder.

`res`: The directory for storing the overlap matrices and NACs.

`min_band`: The minimum KS orbital index to be considered.

`max_band`: The maximum KS orbital index to be considered.

`ks_orbital_homo_index`: The HOMO index for KS orbitals.

`MO_images_directory`: The directory where the molecular orbital isosurfaces images are stored.

`states_to_be_plotted`: The index of the states to be plotted by VMD. If there are many states considered for plotting they should be separated by comma.

`job_init_step`, `nsteps_this_job`, `njob`: Please leave these variables as they are. The program will automatically recognize and fill them.

Now that we have set some of the variables we need to run the python code as `python -c`. First we need to import the function `step2_many_body` from `libra_py.workflows.nbra`. Then we have to set up some variables as a dictionary in `params`. The variables in `params` are as follows:






