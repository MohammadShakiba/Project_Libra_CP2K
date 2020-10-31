# Compute the distributions of the NACs for all states and all steps over the trajectory. Consider the SP and MB bases

Below is an outline of the script nac_dist.py

## 1. Read in the Hvib matrices from the step3/mixed_electron_hole folders.  
## 2. # 2. Divide up into sub-trajectories. These are to be considered our independent nuclear sub-trajectories
In this case, we just take the entire trajectory as a sub-trajectory. This is why we set subtraj_len to = params["nfiles"]

##3.  Obtain the data for the sub trajectory (which in this case is just the data itself) and compute the NAC distributions

Described here is for the case of MB, but the SP case follows analogously. For each pair of electronic states (for all i and js) obtain the nac magnitude in eV. We place this value into a 1x1 CMATRIX object. Please note that one can set the value for the nac cutoff, which is used at 3 different values in our study.
After obtaining all of the nac values, we compute the probability density distributions using data_stat.cmat_distrib()

We then plot the probability density as a function of the bin support using matplotlib

