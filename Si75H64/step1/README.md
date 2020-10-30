# Ab initio Molecular Dynamics with CP2K

The files here contain the inputs required to run the molecular dynamics (MD) in CP2K. This input runs MD for 4000 steps with a time step of 1 fs and a time constant of 16 fs for (CdSe)33 nanocrystal.

`BASIS_MOLOPT` file contains the basis set.
`POTENTIAL` file contains the pseudopotentials for each element.
`dftd3.dat` file contains infomation for considering the van der Waals forces.

Executing CP2K is dependent on the compiled version. You need to load it through use of `module load` or specify the path to the executable CP2K folder:

`export PATH=/full/path/to/executable/cp2k/folder:$PATH`

Here is an example of how to run the MD. The following can be added to your pbs/slurm file for submission:

`mpirun -np 16 cp2k.popt -i cp2k_MD_Si75H64.inp -o output_cp2k_MD_Si75H64.log`
