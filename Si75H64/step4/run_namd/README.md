# Step 4 - Running nonadiabatic molecular dynamics

In this step we run nonadiabatic dynamics (NAD) for the systems in each of the four considered bases. The NAD is run via using the `namd.py` file, which is parallelized for the initial times to increase the speed of calculations. 

## 1. Reading the vibronic Hamiltonian files, the `Hvib` files

We first start by reading the vibronic Hamiltonians for each basis which are stored in the step3 directory. We first set up the required parameters to read such files into a dictionary called `params`. The variables needed for `params` are as follwos:

`data_set_paths`: The path to `res` directory.

`Hvib_re_prefix`, `Hvib_im_prefix`: The prefix used to store the `Hvib`s.

`Hvib_re_suffix`, `Hvib_im_suffix`: The suffix used to store the `Hvib`s.

`nfiles`: Number of files to read.

`nstates`: The number of states to be considered.

`init_times`: The initial time step. Ex) Hvib_sd_10_re has a time step of 10 

`active_space`: The active space which is defined as `list( range(params["nstates"]) )`. The indexing starts from **0**.

Then the `Hvib`s are read using the `step4.get_Hvib2(params)` which is imported from `libra_py.workflows.nbra.step4`.

The same procedure is applied for `Hvib`s of the SD basis. 

## 2. Setting up the initial conditions

We first divide the trajectory up into sub-trajectories. These are to be consdiered our nuclear sub-trajectories. The variables for this part are defined as follows:

`subtraj_time_info`: list of MD trajectory indexes and initial times. Ex) [ [0, 0], [0, 250] ]
                     0th MD trajectory, step 0 and step 250

`nsubtrajs`: The number of sub-trajectories.

`subtraj_len`: The length of each sub-trajectory. Ex) 2000 for [0, 250] would yield the 0th MD trajectory from step 250 to 2250

`params["dt"]`: The time step in the MD in **_atomic_** units. We recommend using the `units.fs2au` which is imported by `from libra_py import units as units` at the beginning of the file.

## 3. Run the dynamics 

Then For each initial condition that we have set up we perform NAD. For each sub-trajectory the input variables are defined as follows:

`params["T"]`: The temperature of the system.
    
`params["ntraj"]`: Number ofsurface hopping stochastic realizations.

`params["sh_method"]`: The surface hopping method.

`params["Boltz_opt"]`: 

`params["nsteps"]`: The number of steps for this dynamics which is set to `subtraj_len`.
 
 `params["init_times"]`: The starting time for the dynamics in the current sub-trajectory.
 
 `params["istates"]`: The initial excited state to start from. This parameter is set with the parameter `istates` list.
 
 `params["decoherence_method"]`: The decoherence method to be considered. 0) FSSH (no decoherence). 1) ID-A. 2) mSDM. 3) DISH 
 
 `params["outfile"]`: The surface hopping output.
 
 `params["nstates"]`: The total number of states to be considered.
 
 We run the NAD through the last lines of the code. Since we have defined the NAD run for each sub-trajectory in the `myfun` function we can run it through use of multiprocessing library of Python with the command `pool.map`. First we create a pool of processors using `pool.map(nprocs)`. Note that the number of processor is optimized to be the same as the number of sub-trajectories. Then, we run the `myfun` function for a set of variables, which here are the subtrajectories and the list of variables becomes `list(range(nsubtrajs))`, using `pool.map( myfunc, list(range(nsubtrajs)) )`. Finally, after running each function is done we close the pool to stop overflow of the computing system using `pool.close()` and `pool.join()`.
 
 
**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs. 
 
 
**_NOTE_:** We highly recommend the user to use the latest inputs used in a recent project about the role of crystal symmetry in nonadiabatic dynamics of CsPbI3 perovskite. The new inputs have better functionality and can be found in [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP).
