# Step 4 - Running nonadiabatic molecular dynamics

In this step we run the nonadiabatic dynamics (NAD) for the systems in the many-body (MB) and single-particle (SP) basis. The input file for NAD is prepared in `namd.py` file which is parallelized for the initial condiotns to increase the speed of calculations. 

## 1. Reading the vibronic Hamiltonian files, the `Hvib` files

We first start by reading the vibronic Hamiltonians which are stored in the `res` directory. We first set up the required parameters into a dictionary called `params`. The variables needed for `params` are as follwos:

`data_set_paths`: The path to `res` directory.

`Hvib_re_prefix`, `Hvib_im_prefix`: The prefix used to store the `Hvib`s.

`Hvib_re_suffix`, `Hvib_im_suffix`: The suffix used to store the `Hvib`s.

`nfiles`: Number of files to read.

`nstates`: The number of states to be considered.

`init_times`: The initial time step to start reading the files.

`active_space`: The active space which is defined as `list( range(params["nstates"]) )`. The indexing starts from **0**.

Then the `Hvib`s are read using the `step4.get_Hvib2(params)` which is imported from `libra_py.workflows.nbra.step4`.

The same procedure is applied for `Hvib`s of the SD basis. 

## 2. Setting up the initial conditions

We first divide the trajectory up into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories. The variables for this part are defined as follows:

`nsubtrajs`: The number of sub-trajectories.

`subtraj_len`: The length of each sub-trajectory.

`start_time`: The start time for the first sub-trajectory.

`subtraj_increment`: The subtrajectories increment size.

`params["dt"]`: The time step in the MD in **_atomic_** units. We recommend using the `units.fs2au` which is imported by `from libra_py import units as units` at the beginning of the file.

## 3. Run the dynamics 

Then For each initial condition that we have set up we perform NAD. For each sub-trajectory the input variables are defined as follows:

`params["T"]`: The temperature of the system.
    
`params["ntraj"]`: Number ofsurface hopping stochastic realizations.

`params["sh_method"]`: The surface hopping method.

`params["Boltz_opt"]`: 

`params["nsteps"]`: The number of steps for this dynamics which is set to `subtraj_len`.
 
 `params["init_times"]`: The starting time for the dynamics in the current sub-trajectory.
 
 `params["istate"]`: The initial excited state to start from. This parameter is set with the parameter `istates` list.
 
 `params["decoherence_method"]`: The decoherence method to be considered. 
 
 `params["outfile"]`: The surface hopping output.
 
 `params["nstates"]`: The total number of states to be considered.
 
 Then, the NAD is run through `step4.run(hvib, params)`. Here we use the `step4.run( [ hvib_mb_subtrajs[ subtraj ] ], params)` for each sub-trajectory. The first parameter of the `step4.run` are the `Hvib`s which should be represented as a list of lists.
 
 We run the NAD through the last lines of the code. Since we have defined the NAD run for each sub-trajectory in the `myfun` function we can run it through use of multiprocessing library of Python with the command `pool.map`. First we create a pool of processors using `pool.map(nprocs)`. Note that the number of processor is optimized to be the same as the number of sub-trajectories. Then, we run the `myfun` function for a set of variables, which here are the subtrajectories and the list of variables becomes `list(range(nsubtrajs))`, using `pool.map( myfunc, list(range(nsubtrajs)) )`. Finally, after running each function is done we close the pool to stop overflow of the computing system using `pool.close()` and `pool.join()`.
 
 
 
 
 

