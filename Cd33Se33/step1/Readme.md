# Ab initio Molecular Dynamics with CP2K

The files herre contain the inputs required to run the molecular dynamics (MD) in CP2K. This input runs MD for 4000 steps with a time step of 1 fs and a time constant of 16 fs.

`BASIS_MOLOPT` file contains the basis set.
`POTENTIAL` file contains the pseudopotentials for each element.
`dftd3.dat` file contains infomation for considering the van der Waals forces.

To run the MD you can add the following command to your pbs/slurm file:

`mpirun -np 16 cp2k.popt -i cp2k_MD_Cd33Se33.inp -o output_cp2k_MD_Cd33Se33.log`
