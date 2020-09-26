# Step 3 - Forming the Slater determinant basis NACs from KS overlaps

In this step we form the vibronic Hamiltonian matrices of the Slater determinant (SD) basis. The steps for doing this task are as follows:

## 1. Read energies and time-verlaps of Kohn-Sham basis

We first start by reading the files that have the energies and time-overlap matricies in the Kohn-Sham (KS) basis i.e. `E_ks` and `St_ks`. We start by setting the `params` dictionary. First we set up the variable `data_set_paths` which is the path to `St_ks` and `E_ks` files which are located in the `res` directory. Other variables are as follows:

`data_dim`: The number of KS states which is the number of rows in the `E_ks` files.

`active_space`: The active space for KS orbitals. Here, we set it to `range(data_dim)`.

`start_time`: The initial time step.

`finish_time`: The final time step plus one (`+1`).

After setting up the `params` variable we start reading the data for both KS energies and overlaps using `data_read.get_data_sets(params)`.


## 2. Obtain the unique SDs from all TD-DFT calculations over the whole trajectories

The unique SD are obtained from running `get_unique.py`. This file will read the slurm (`pbs` or `bash`) output file (that is, output by Libra) that was made from processing the output files of the TD-DFT CP2K calculations and extracts the unique set of SDs over all timesteps. These SDs are given in terms of their 1-electron KS excitation from a homo index w.r.t **1**. We then **manually** add the ground state to the `unique_SDs`. For example `['28, 28']` as is shown in the `step3.py` file, where 28 is the `homo index` from **1**. This show some excitaiton from 28 to itself, which results in no change, and is thus the HOMO. 

Then, we transform the unique SD excitatons as given by CP2K into the format that Libra expects. We considered only transitions within the spin-restricted reference, so our excitations are of 1 spin type. Here we choose the alpha electrons only. 

## 3. Oredering the SDs by their energies

Now, for each timestep, we compute the energies of the SDs and sort them by energy. This prevents the state energies from crossing in this basis. However, note that the St_ks files (time-overlaps in the KS basis) have no knowledge of the reordering of the SDs by energies.

## 4. Computing the mid-point energies and time-overlaps according to the energy ordered SDs

Note that by computing wavefunciton overlaps of the energy ordered SDs using the properties of the KS orbtials results in the time-overlaps in the SD basis taking their largest magnitudes in off-diagonal locations. This is because the ordering of the wavefuncitons in the KS basis is not consistent with the ordering of the SD basis based on energy ordering. Only in special cases does the energy ordering of the SD basis match the ordering of the KS basis (such as for `homo->lumo+N` (electron-only) or `homo-N->lumo` (hole-only) excitations). Therefore, we will need to apply a state reordering procedure to the `St_sd` matricies computed here.

## 5. Applying the state reordering procedured 

In the nes step we apply the state reordering as was mentioned above, followed by a correction to the phases.

## 6. Form `Hvib_sd`

Now we form the vibronic Hamiltonian from the `St_sd` using the Hammes-Schiffer-Tully method.

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs.
