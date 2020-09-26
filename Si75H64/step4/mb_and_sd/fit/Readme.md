# Step 4 - Computing the timescales by fitting the NAD results by an stretched-compressed exponential function:

In order to fit the NAD results to an stretched-compressed exponential function we first need to define it. There are various function available but we have adjusted the working fitting code to the stretched-compressed exponential one which is defined as in `func3` in the code.

We first set up the paramters for reading/sorting the data. The time step is defined as `dt` in femtosecond. Other variables are as follows:

`basis_options`: This list includes the dynamic basis. Here we have `sd` and `mb`. In fact, these are defined in the `namd.py` when outputting the NAD results.

`decoherence_options`: This list includes the decoherence options which were also done in the `namd.py`. The list input will be used to identify the NAD output files. For example we have used the following decoherence option in the `namd.py`: `fssh`, `ida`, `msdm`

`decoherence_options_names`: The names of the decoherence options used for plotting (the names will be set as titles).

`dynamics_option`: This variable defines for the fit function to be applied for which state population. There are different options here: **`0`** yields the hot energy decay timescales, **`1`** yields the recovery timescales, and **`2`** is for recombination timescales. For the recovery timescales, we need to define the recovery of which states to be considered. Here, we consider the recovery of the lowest excited states for CdSe NC. The user can define lowest two or three excited states by subtracting the population of those states as well. For example, the `y = 1.0 - ( float( namd_data_line[ 6 ] ) )` in the code defines the recovery of the first excited states. One can obtain the recovery of the first two excited states by `y = 1.0 - ( float( namd_data_line[ 6 ] ) - float( namd_data_line[ 9 ] ) )` in which `float( namd_data_line[ 6 ] )` shows the population of the first excited state and `float( namd_data_line[ 9 ] )` shows the population of the second excited state in the NAD out files.

`nsubtrajs`: The number of sub-trajectories. This is needed because we want to obtain the total average of the timescales as well.

`istates`: The initial excited states as in the `namd.py` file. This input is also required to identify the output files by NAD calculations.

`r_squared_tol`: Here we only consider the fits that have an R squared factor of larger than a specific value. This value is defined by `r_squared_tol`.

After fitting, we obtain the average and the standard deviation value of all fitted data. Then, we use the margin of error for a sample proportion as is shown [here](https://barbatti.org/2018/04/18/how-many-trajectories-should-i-run/).


The detailed results are output by `python fit.py`. One can obtain the average values and their corresponding error bars from the output using the commad `python fit.py | grep "Reading\|final"`. Other alternatives would be to do first `python fit.py > output.log` and then using `grep` to find the needed values by `cat output.log | grep -i "Reading\|final\|beta"` which `-i` is for case insensitive while using `grep`.

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs.
