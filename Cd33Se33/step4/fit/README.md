# Step 4 - Computing the timescales by fitting the NAD results via a stretched-compressed exponential function

The python script used to perform the fitting is the file "fit_template.py"

## First, we define two funcitons that will be useful during the fitting process
The first two functions defined in this file are the fitting function "func" and a wrapper function to the actual fitting. The wrapper function "fit_data" returns the fit of the data as well as the optimal fit parameters and the r squared values of the fits

## Second, we set up the paramters for reading/sorting the data. The time step is defined as `dt` in femtosecond. Other variables are as follows:

`decoherence_options`: This list includes the decoherence options which were also done in the `namd.py`. The list input will be used to identify the NAD output files. For example we have used the following decoherence option in the `namd.py`: `fssh`, `ida`, `msdm`

`decoherence_options_names`: The names of the decoherence options used for plotting (the names will be set as titles).

'B_max': This is a global variable simply set to a very low number. In the stretch-compress fit, one can optionall include and shift tot he function. Here, we set it to a very low value to essentially not include it.

`dynamics_option`: This variable defines the type of dynamics process we want to consider. There are different options here: **`0`** yields the hot energy decay timescales, **`1`** yields the recovery timescales of some user defined electronic states, and **`2`** is for recombination process. For the recovery timescales, we need to define to which states the recovery process is to be considered. 

`nsubtrajs`: The number of sub-trajectories. This is needed because we want to obtain the total average of the timescales as well. This number corresponds to the number of initial times considered in the NAMD

`taus, betas, dyns`: These are empty lists which will be populated with data from the fittings

`basis_options`: which basis do we want to consider? mb = many-body, etc. Please keep them in this order 

`r_squared_tol`: Here we only consider the fits that have an R squared factor of larger than a specific value. This value is defined by `r_squared_tol`.

`istates`: The initial excited states as defined in the `namd.py` file.

After fitting, we obtain the average and the standard deviation value of all fitted data. Then, we use the margin of error for a sample proportion as is shown [here](https://barbatti.org/2018/04/18/how-many-trajectories-should-i-run/).


The detailed results are output by `python fit.py`. One can obtain the average values and their corresponding error bars from the output using the commad `python fit.py | grep "Reading\|final"`. Other alternatives would be to do first `python fit.py > output.log` and then using `grep` to find the needed values by `cat output.log | grep -i "Reading\|final\|beta"` which `-i` is for case insensitive while using `grep`.

**_NOTE_:** - Please note that the paths currently defined in these files may not be the correct paths for you. Please adjust all paths to your specific needs.


**_NOTE_:** We highly recommend the user to use the latest inputs used in a recent project about the role of crystal symmetry in nonadiabatic dynamics of CsPbI3 perovskite. The new inputs have better functionality and can be found in [this link](https://github.com/AkimovLab/Project_CsPbI3_MB_vs_SP).
