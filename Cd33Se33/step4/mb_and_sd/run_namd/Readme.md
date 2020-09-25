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

`params["T"]`: Temperature                  = 300.0
    params["ntraj"]              = 250
    params["sh_method"]          = 1
    params["Boltz_opt"]          = 1

    params["nsteps"] = subtraj_len
    params["init_times"] = [0]

