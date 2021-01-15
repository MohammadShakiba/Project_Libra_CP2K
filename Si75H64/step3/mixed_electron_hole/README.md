# Step 3 - Forming the vibronic Hamiltonian in the many-body and Slater determinant basis

In this step we form the vibronic Hamiltonian in the many-body and Slater determinant (SD) basis. The steps for doing this task are as follows:

## 1. Read energies and overlaps computed in the Kohn-Sham basis

We first start by reading the files that have the energies and time-overlap matricies in the Kohn-Sham (KS) basis i.e. `E_ks`, `S_ks` and `St_ks`. We start by setting the `params` dictionary. First we set up the variable `data_set_paths` which is the path to `S_ks`, `St_ks`, and `E_ks` files, which are located in the `res` directory. Other variables are as follows:

`data_dim`: The number of KS states which is the number of rows in the `E_ks` files.

`active_space`: The active space for KS orbitals. Here, we set it to `range(data_dim)`.

`start_time`: The initial time step.

`finish_time`: The final time step plus one (`+1`).

After setting up the `params` variable we start reading the data for both KS energies and overlaps using `data_read.get_data_sets(params)`. We recommend to set finish_time to 2 and to uncomment the sys.exit(0) after the St_ks section in order to test the reading of the KS data. If matricies containing all zeros prints to the screen, the reading of the data was not done correctly. 

## 2. Obtain the excitation energies and ci coefficients or all timesteps
Make a directory called "all_logfiles" and then run the following command: 
`for file in $(find ../../step2/wd -name '*.log'); do cp $file all_logfiles/.; echo $file; done`
Please note, that you must provide the path to the wd directory in the step2 calculations after the `find` command

## 3. Process the unique SD bases, which involves the following
3.1. Reindexing the SD bases from the native ES format to what Libra expects
3.2. Sorting the SD bases at each timestep by their energies 

Note: The way in which SD states are to be written has been reformulated in the most recent versions of Libra. To reproduce the results of this study, one needs to set
the input parameter do_sort=True in the funciton sd2indx in the file $path/libra_py/workflows/nbra/mapping.py and express the SDs in the following form:
[ all alpha orbtials, all beta orbitals ], where the alpha and beta orbtials correspond to row or column elements (from 1) of the spin-block matrix in the KS absis. Beta electrons shall have a
negative sign in fron of them. Ex) For 10 alpha and 10 beta orbitals, with HOMO = 5 (from 1), the GS will be [ 1, 2, 3, 4, 5, -10, -11, -12, -13, -14 ]

## 4. We need to make the linear transformation matrices. This involes:
4.1. Make the list of ci_coefficients for each step in the way Libra accepts
4.2. Making the T matrix from this list keeping in mind the SD bases

## 5. Use the SD to CI transformation matrices to convert from St_sd -> St_ci
We also need to make the CI energy matrix from the excitation energies and finally make the Hvib in the CI basis
the SD basis (Hvib, overlaps) is output at the end as well     

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs.
Also, this script is quite "raw" and is a result of the "brand-new: nature of the NAMD in MB basis feature. In subsuent versions of Libra, it is possible than many aspect of this script will be refined
and placed into either libra_py or the Libra's core modules. Also, due to the large number of files processed and potentially large number of electrons in the SDs, it is possible this script may take several hours. It is recommend to submit this file for computation. The file is currently parallelized with 24 processors via `mp.pool()`. 


**_NOTE_:** We highly recommend the user to use the latest inputs used in a recent project about the role of crystal symmetry in nonadiabatic dynamics of CsPbI3 perovskite. The new inputs have better functionality and can be found in [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP).
