from libra_py import data_stat
from libra_py import CP2K_methods
import os
import sys
import multiprocessing as mp
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.colors
import math
import time
#plt.switch_backend('agg')

# from libra_py import CP2K_methods
params = { }
params["number_of_states"] = 30
params["tolerance"] = 0.01
params["isUKS"] = 0
ci_coeffs = []

logfiles = glob.glob('all_logfiles/*.log')
#logfiles = glob.glob('optimized_logfile/*.log')

for logfile in logfiles:
    params.update({"logfile_name": logfile})
    excitation_energies, ci_basis, ci_coefficients, spin_components = CP2K_methods.read_cp2k_tddfpt_log_file(params)
    for j in range(len(ci_coefficients)):
        for k in range(len(ci_coefficients[j])):
            ci_coefficients[j][k] = ci_coefficients[j][k]**2
    ci_coeffs.append(ci_coefficients)

print(ci_coeffs)
print(len(ci_coeffs))
print(len(ci_coeffs[0]))

# The goal is to get a list of len nstates. Each element is a list of lists of size nsteps
nsteps = len(ci_coeffs)
nstates = params["number_of_states"]
nsds = 3
coeffs = []
coeffs_avg   = []
coeffs_error = []

plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
plt.subplot(1,1,1)

plt.title("300 K (CdSe)$_{33}$", fontsize=10)
#plt.title("0 K (CdSe)$_{33}$", fontsize=10)

plt.xlabel('State Index', fontsize=10)
plt.ylabel('c$_{i}^2$',   fontsize=10)
plt.xticks([0,10,20,30])

for state in range(nstates):

    coeffs.append( [] )
    coeffs_avg.append( [] )
    coeffs_error.append( [] )

    for sd in range( nsds ):

        coeffs[state].append( [] )
        coeffs_avg[state].append( [] )
        coeffs_error[state].append( [] )

        for step in range( nsteps ):
            if len( ci_coeffs[step][state] ) < nsds and sd > len( ci_coeffs[step][state] )-1:
                coeffs[state][sd].append( 0.0 )
            else:
                coeffs[state][sd].append( ci_coeffs[step][state][sd] )
     
        mb_coeff_avg, mb_coeff_std = data_stat.scalar_stat( coeffs[state][sd] )
        coeffs_avg[state][sd].append( mb_coeff_avg )
        coeffs_error[state][sd].append( 1.96 * mb_coeff_std / math.sqrt(nsteps) )

        if sd == 0:
            print("std = ", mb_coeff_std)     

        if sd == 0:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='green', markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        elif sd == 1:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='blue',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')

        else:
            plt.plot(     state+1, mb_coeff_avg, color="black", marker='s', markerfacecolor='red',  markeredgewidth=0.4, markersize=5)
            plt.errorbar( state+1, mb_coeff_avg, yerr=coeffs_error[state][sd], linestyle="None", color='black')
 

plt.tight_layout()
plt.savefig('Cd33Se33_excitation_analysis_300K.png', dpi=300)
#plt.savefig('Cd33Se33_excitation_analysis_0K.png', dpi=300)

