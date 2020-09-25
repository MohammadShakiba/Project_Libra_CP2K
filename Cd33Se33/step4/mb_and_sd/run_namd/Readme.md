# Step 4 - Running nonadiabatic molecular dynamics

In this step we run the nonadiabatic dynamics (NAD) for the systems in the many-body (MB) and single-particle (SP) basis. The input file for NAD is prepared in `namd.py` file which is parallelized for the initial condiotns to increase the speed of calculations. 

## 1. Reading the vibronic Hamiltonian files, the `Hvib`s

We first start by reading the vibronic Hamiltonians which are stored in the `res` directory. We first set up the required parameters into a dictionary called `params`. The variables needed for `params` are as follwos:

`data_set_paths`: The path to `res` directory.

`Hvib_re_prefix`, `Hvib_im_prefix`: The prefix used to store the `Hvib`s.

`Hvib_re_suffix`, `Hvib_im_suffix`: The suffix used to store the `Hvib`s.

`nfiles`: Number of files to read.

`nstates`: The number of states to be considered.

`init_times`: The initial time step to start reading the files.

`active_space`: The active space which is defined as `list( range(params["nstates"]) )`. The indexing starts from **0**.

Then the `Hvib`s are read using the `step4.get_Hvib2(params)` which is imported from `libra_py.workflows.nbra.step4`.


