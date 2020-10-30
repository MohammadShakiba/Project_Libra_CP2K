# Step 3 - Forming the Slater determinant basis NACs from KS overlaps

In this step we form the vibronic Hamiltonian matrices of the Slater determinant (SD) basis. The steps for doing this task are as follows:

## 1. Read energies and time-verlaps of Kohn-Sham basis

We first start by reading the files that have the energies and time-overlap matricies in the Kohn-Sham (KS) basis i.e. `E_ks` and `St_ks`. We start by setting the `params` dictionary. First we set up the variable `data_set_paths` which is the path to `St_ks` and `E_ks` files which are located in the `res` directory. Other variables are as follows:

`data_dim`: The number of KS states which is the number of rows in the `E_ks` files.

`active_space`: The active space for KS orbitals. Here, we set it to `range(data_dim)`.

`start_time`: The initial time step.

`finish_time`: The final time step plus one (`+1`).

After setting up the `params` variable we start reading the data for both KS energies and overlaps using `data_read.get_data_sets(params)`.

**_NOTE_:** - you must make a directory called "all_logfiles" and then run the following command 
for file in $(find ../step2/wd -name '*.log'); do cp $file all_logfiles/.; echo $file; done

Now, for each timestep, the energies of the SDs is computed and sorted by their energies. This prevents the state energies from crossing in this basis. However, note that the St_ks files (time-overlaps in the KS basis) have no knowledge of the reordering of the SDs by energies.

## 4. Computing the mid-point energies and time-overlaps according to the energy ordered SDs

Note that by computing wavefunciton overlaps of the energy ordered SDs using the properties of the KS orbtials results in the time-overlaps in the SD basis taking their largest magnitudes in off-diagonal locations. This is because the ordering of the wavefuncitons in the KS basis is not consistent with the ordering of the SD basis based on energy ordering. Only in special cases does the energy ordering of the SD basis match the ordering of the KS basis (such as for `homo->lumo+N` (electron-only) or `homo-N->lumo` (hole-only) excitations). Therefore, we will need to apply a state reordering procedure to the `St_sd` matricies computed here.

## 5. Form `Hvib_sd`

Now we form the vibronic Hamiltonian from the `St_sd`

One then computes the transformation matrix "SD2CI" which converts the SD basis to the CI-like (MB) basis. After transformation to the MB basis, phase corrections and state-reordering are applied. Next, the Hamiltonian
in the MB basis is made

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs.
