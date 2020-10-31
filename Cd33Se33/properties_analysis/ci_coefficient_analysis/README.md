# Plot the average CI coefficients over the molecular dynamics trajectory

Here we plot the average CI coefficients obtained from TD-DFPT calculations for each excited states. To this end we need all the log files obtained from CP2K TD-DFPT calculations. So, we need
to copy all the log files in one folder you can use the following bash script:

```
mkdir all_logfiles
for file in $(find path/to/wd -name '*.log'); do cp $file all_logfiles/.; echo $file; done
```

Then we need to specify the parameters as follows:

`params["number_of_states"]`: The number of excited states considered in the TD-DFPT calculations. Here we use _30 _ for CdSe NC.

`params["tolerance"]`: The tolerance factor for CI coefficients.

`params["isUKS"]`: The spin-polarization flag. For spin restricted case we use _0_.

The we need to find and append the names of all log files into one variable so that we can use the multiprocessing feature which is appended in the variable `logfiles`, using 
`logfiles = glob.glob('all_logfiles/*.log')`. After that, we can take the average of the CI coefficients and plot them.
