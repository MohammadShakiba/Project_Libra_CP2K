# Step 3 - Forming the Slater determinant basis from KS overlaps





## 1. Read the files that have the energies and time-overlap matricies in the Kohn-Sham basis, E_ks and St_ks.

## 2. shown below are excitations. These are obtained from running "get_unique.py", which goes into the slurm output file (that is, output by Libra) that was made from processing the output files of the TD-DFT cp2k calculations and extracts the unique set of SDs over all timesteps. These SDs are given in terms of their 1-electron KS excitation from a homo index w.r.t 1. We then manually added the ground state as ['28, 28'], where 28 is the homo index from 1. This show some excitaiton from 28 to itself, which results in no change, and is thus the HOMO. See section 2.1.

## 3. Now, for each timestep, we compute the energies of the SDs and sort them by energy. This prevents the state energies from crossing in this basis. However, note that the St_ks files (time-overlaps in the KS basis) have no knowledge of the reordering of the SDs by energies

## 4. Compute the mid-point energies and time-overlaps according to the energy ordered SDs. Note that by computing wavefunciton overlaps of the energy ordered SDs using the properties of the KS orbtials results in the time-overlaps in the SD basis taking their largest magnitudes in off-diagonal locations. This is because the ordering of the wavefuncitons in the KS basis is not consistent with the ordering of the SD basis based on energy ordering. Only in special cases does the energy ordering of the SD basis match the ordering of the KS basis (such as for homo->lumo+N or homo-N->lumo only excitations). Therefore, we will need to apply a state reordering procedure to the St_sd matricies computed here.

## 5. Now, apply the state reordering procedured mentioned above, followed by a correction to the phases.

## 6. Form Hvib_sd
