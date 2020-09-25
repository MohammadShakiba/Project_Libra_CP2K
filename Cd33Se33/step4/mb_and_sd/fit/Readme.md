# Step 4 - Computing the timescales by fitting the NAD results by an stretched-compressed exponential function:

In order to fit the NAD results to an stretched-compressed exponential function we first need to define it. There are various function available but we have adjusted the working fitting code to the stretched-compressed exponential one which is defined as in `func3` in the code.

We first set up the paramters for reading/sorting the data. The time step is defined as `dt` in femtosecond. Other variables are as follows:

`basis_options`: This list includes the dynamic basis. Here we have `sd` and `mb`. In fact, these are defined in the `namd.py` when outputting the NAD results.

`decoherence_options`: This list includes the decoherence options which were also done in the `namd.py`. The list input will be used to  ["fssh","ida", "msdm"]#, "dish"]

decoherence_options_names = ["FSSH","ID-A","mSDM"]

fit_options = [3]

# 0 = "Hot_energy decay"
# 1 = "Recovery"
# 2 = "Recombination"
dynamics_option = 1

nsubtrajs = 21
# lower: sd, mb
#istates = [ [5,6,7,8], [5,6,7,8] ]
# upper: sd, mb
istates = [ [14,15,16,17,18], [14,15,16,17,18] ]

r_squared_tol = 0.8
