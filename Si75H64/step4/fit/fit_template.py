#***********************************************************
# * Copyright (C) 2020 Brendan Smith, Mohammad Shakiba, Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/
import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import data_stat
from libra_py import units
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.optimize import curve_fit





####################
# 1. Define the needed helper functions. This includes the fitting function and a wrapper funciton to call it
def func(t, tau, beta, E0, B):
    return E0*np.exp( -(t/tau)**beta ) + B

def fit_data( xdata, ydata, B_max ):

    ydata_fit = []
    popt, pcov = curve_fit( func, xdata, ydata, bounds=([0.0, 0.0, ydata[0]-0.001, 0.0], [np.inf, np.inf, ydata[0]+0.001, B_max]))
    residuals  = ydata - func(xdata, *popt)
    ss_res, ss_tot = np.sum(residuals**2), np.sum((ydata - np.mean(ydata))**2)
    r_squared  = 1.0 - (ss_res / ss_tot)
    tau, beta, E0, B = popt

    # If we have no population recovery, our ydata will be a list of 1.0 values. 
    # However, the plotting funciton will error and can even give something looking like decay.
    # So, in this case, just return a list of 1.0s
    if all(ydata_el==1.0 for ydata_el in ydata) == True or r_squared < 0:
        for i in range(len(xdata)):
            ydata_fit.append(1.0)
    else:
        for i in range(len(xdata)):
            ydata_fit.append(func(xdata[i], *popt))

    # Uncomment the below section to plot the fit. This can be done to check a special case of intrest, if need be 
    """
    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title("TEST", fontsize=8) 
    plt.xlabel('Time, fs',   fontsize=8)
    plt.plot(xdata, ydata,     linewidth=1.5, color="black", label="data")
    plt.plot(xdata, ydata_fit, linewidth=1,   color="green", label="fit") 
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("tmp.png")
    """

    return ydata_fit, tau, beta, E0, B, r_squared



####################
# 2. Set paramters for reading / sorting the data
dt = 0.5 #fs
decoherence_options       = ["fssh","ida", "msdm"]
decoherence_options_names = ["FSSH","ID-A","mSDM"]
fit_options = [3]
B_max = 0.001
# 0 = "Hot_energy decay"
# 1 = "Recovery"
# 2 = "Recombination"
dynamics_option = 0
nsubtrajs = 8
taus  = []
betas = []
dyns  = []
decoherence_option_count = 0
# Hot energy
if dynamics_option == 0:
    basis_options = ["mb","mixed_sd","elec_sd","hole_sd"]
    r_squared_tol = 0.8
    istates = [ [45,46,47,48], [45,46,47,48], [7,8,9], [13,14,15] ]
# Recovery
if dynamics_option == 1:
    basis_options = ["mb","mixed_sd","elec_sd","hole_sd"]
    r_squared_tol = 0.8
    istates = [ [45,46,47,48], [45,46,47,48], [7,8,9], [13,14,15] ]




#####################
# 3. Enter the "ultimate" fitting loop. For each decoherence option, we will perofrm dynamics using both the SD and MB bases
for decoherence_option in decoherence_options:

    print("Number of initial conditions before R2 filtering", nsubtrajs*len(istates[0]))
    taus.append(  [] )
    betas.append( [] )
    dyns.append(  [] )

    # 3.1. For each basis option choice (SD and MB) we have a number of initial conditions, we will obtain a tau for each
    basis_option_count = 0
    for basis_option in basis_options:
        
        taus[decoherence_option_count].append(  [] )
        betas[decoherence_option_count].append( [] )
        dyns[decoherence_option_count].append(  [] )     

        # 3.2. Now, the initial conditions are for 4 nuclear sub-trajectories, each having 3 initial "istate" conditions
        for subtraj in range(nsubtrajs):

            # 3.3.. Again, each istate is also an initial condition, here. So for each initial condition, we collect a tau
            for istate in istates[basis_option_count]:

                xdata     = []
                ydata     = []
                ydata_fit = []

                filename = "../step4_surface_hopping/_out_si75h64_"+basis_option+"_nbra_namd_"+decoherence_option+"_"+str(subtraj)+"_"+str(istate)+".txt"
                outname  = "_fit_"+basis_option+"_"+decoherence_option+"_subtraj_"+str(subtraj)+"_istate"+str(istate)+".png"

                f = open(filename)
                A = f.readlines()
                sz = len(A)
                f.close()
                count = 0
                #sys.exit(0)

                print ("\nReading file",filename)

                for i in range(sz):

                    namd_data_line = A[i].strip().split()

                    # Hot Energy Decay
                    if dynamics_option == 0:
                        # mb dynamics
                        if basis_option_count == 0:
                            y = float( namd_data_line[ 155 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.00

                        # mixed_sd dynamics
                        elif basis_option_count == 1:
                            y = float( namd_data_line[ 257 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.00

                        # elec dynamics
                        elif basis_option_count == 2:
                            y = float( namd_data_line[ 38 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.00

                        # hole dynamics
                        elif basis_option_count == 3:
                            y = float( namd_data_line[ 68 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.00

                        y = y * units.au2ev
                        #print(y)

                    # Recovery
                    elif dynamics_option == 1:

                        #"""
                        # mb, mixed_sd
                        if basis_option_count == 0 or basis_option_count == 1:
                            population_recovered = float( namd_data_line[ 3 ] ) + float( namd_data_line[ 6 ] ) + float( namd_data_line[ 9 ] )

                        # elec_sd
                        elif basis_option_count == 2:
                            population_recovered = float( namd_data_line[ 3 ] ) + float( namd_data_line[ 6 ] ) + float( namd_data_line[ 9 ] )  

                        # hole_sd
                        elif basis_option_count == 3:
                            population_recovered = float( namd_data_line[ 3 ] ) + float( namd_data_line[ 6 ] ) + float( namd_data_line[ 9 ] ) 
                        #"""

                        #population_recovered = float( namd_data_line[ 3 ] ) + float( namd_data_line[ 6 ] ) 
                        y = 1.0 - population_recovered
                        

                    # Ground state recovery
                    elif dynamics_option == 2:
                            y = 1.0 - ( float( namd_data_line[ 3 ] ) ) 

                    xdata.append(i*dt)
                    ydata.append(y)

                # Accept all subtrajectories for plotting, but only compute the statistics of the taus from 
                # subtrajectories that have r_squared > r-squared_tol
                dyns[decoherence_option_count][basis_option_count].append( np.array( ydata ) )
                try:
                    ydata_fit, tau, beta, E0, B, r_squared = fit_data( xdata, ydata, B_max ) 
                    # Append a tau for each initial condition, here, initial conditions are istates ..
                    if r_squared >= r_squared_tol:
                        taus[decoherence_option_count][basis_option_count].append(tau)
                        betas[decoherence_option_count][basis_option_count].append(beta)
                        print ("r_squared = ", r_squared)
                        print ("Printing optimal fitting paramters:")
                        print ("B    = ", B)
                        print ("E0   = ", E0, "eV")
                        print ("beta = ", beta)
                        print ("tau  = ", tau, "fs")
 
                except:
                    pass              


        niconds = len(taus[decoherence_option_count][basis_option_count])
        print("final ", decoherence_options_names[decoherence_option_count], basis_options[basis_option_count], "niconds", niconds)

        # Sometimes, we do not have any population decay, or any successful fits. in this case, just explain this
        if len( taus[decoherence_option_count][basis_option_count] ) < 1:
            print("No taus or betas were found for ", decoherence_options[decoherence_option_count], basis_options[basis_option_count])
            pass        
        else:
            taus_avg, taus_std = data_stat.scalar_stat( taus[decoherence_option_count][basis_option_count] )       
            error_95_taus = 1.96 * taus_std / math.sqrt(niconds)

            betas_avg, betas_std = data_stat.scalar_stat( betas[decoherence_option_count][basis_option_count] )
            error_95_betas = 1.96 * betas_std / math.sqrt(niconds)
            print ("\nFinished computing data for ", decoherence_options[decoherence_option_count], [basis_option_count])
            print ("Average decay lifetime over every NAD run: ",             taus_avg, " fs")
            print ("Std_dev decay lifetime over every NAD run: ",             taus_std, " fs")
            print ("Average population decay lifetime with error (95% ci): ", taus_avg, "+-", error_95_taus, " fs\n")
            print ("final taus ", decoherence_options[decoherence_option_count], basis_options[basis_option_count],  taus_avg, "+-", error_95_taus,  "fs\n") 
            print ("final betas", decoherence_options[decoherence_option_count], basis_options[basis_option_count], betas_avg, "+-", error_95_betas, "fs\n")

        basis_option_count += 1

    #sys.exit(0)

    ####################
    # 4.0 We have compute the dynamics for all bases (mb, mixed_sd, elec_sd, hole_sd) for a given decoherence scheme. Now, let's make our plots
    ydata_avg_mb = sum(dyns[decoherence_option_count][0]) / len(dyns[decoherence_option_count][0])
    ydata_avg_mb_fit = []

    ydata_avg_mixed_sd = sum(dyns[decoherence_option_count][1]) / len(dyns[decoherence_option_count][1])
    ydata_avg_mixed_sd_fit = []

    ydata_avg_elec_sd = sum(dyns[decoherence_option_count][2]) / len(dyns[decoherence_option_count][2])
    ydata_avg_elec_sd_fit = []

    ydata_avg_hole_sd = sum(dyns[decoherence_option_count][3]) / len(dyns[decoherence_option_count][3])
    ydata_avg_hole_sd_fit = []

    ydata_avg_mb_fit, tau_avg_mb, beta, E0, B, r_squared = fit_data( xdata, ydata_avg_mb, B_max )
    ydata_avg_mixed_sd_fit, tau_avg_mixed_sd, beta, E0, B, r_squared = fit_data( xdata, ydata_avg_mixed_sd, B_max )  
    ydata_avg_elec_sd_fit, tau_avg_elec_sd, beta, E0, B, r_squared   = fit_data( xdata, ydata_avg_elec_sd, B_max )
    ydata_avg_hole_sd_fit, tau_avg_hole_sd, beta, E0, B, r_squared   = fit_data( xdata, ydata_avg_hole_sd, B_max )    

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title("Si$_{75}$H$_{64}$ "+decoherence_options_names[decoherence_option_count], fontsize=12)
    plt.xlabel('Time, fs',        fontsize=12)
    plt.xticks(fontsize=12) 
    if dynamics_option == 0:
        plt.ylabel('Hot Energy, eV', fontsize=12)
        plt.yticks([0.0,0.3,0.6,1.0],fontsize=12)
        plt.ylim(0,1.0)
    elif dynamics_option == 1:
        plt.ylabel('Population',     fontsize=12)
        plt.yticks([0.0,0.3,0.6,1.0],fontsize=12)
        plt.ylim(0,1.0)
    for i in range( len(dyns[decoherence_option_count][1]) ):
        plt.plot( xdata, dyns[decoherence_option_count][1][i], label="", linewidth=1, color = "pink")
    for i in range( len(dyns[decoherence_option_count][0]) ):
        plt.plot( xdata, dyns[decoherence_option_count][0][i], label="", linewidth=1, color = "skyblue")
    plt.plot(xdata, ydata_avg_mixed_sd_fit, linewidth=3, color="red",   label="Single-particle")
    plt.plot(xdata, ydata_avg_mb_fit,       linewidth=3, color="blue",      label="Many-body")
    plt.tight_layout()

    if dynamics_option == 0:
        plt.savefig("Si75H64_mb_mixed_sd_Hot_Energy_Decay_"+decoherence_options[decoherence_option_count])
    elif dynamics_option == 1:
        plt.savefig("Si75H64_mb_mixed_sd_Recovery_"+decoherence_options[decoherence_option_count])

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title("Si$_{75}$H$_{64}$ "+decoherence_options_names[decoherence_option_count], fontsize=12)
    plt.xlabel('Time, fs',        fontsize=12)
    plt.xticks(fontsize=12)
    if dynamics_option == 0:
        plt.ylabel('Hot Energy, eV', fontsize=12)
        plt.yticks([0.0,0.3,0.6,1.0],fontsize=12)
        plt.ylim(0,1.0)
    elif dynamics_option == 1:
        plt.ylabel('Population',     fontsize=12)
        plt.yticks([0.0,0.3,0.6,1.0],fontsize=12)
        plt.ylim(0,1.0)
    for i in range( len(dyns[decoherence_option_count][3]) ):
        plt.plot( xdata, dyns[decoherence_option_count][3][i], label="", linewidth=1, color = "mediumorchid")
    for i in range( len(dyns[decoherence_option_count][2]) ):
        plt.plot( xdata, dyns[decoherence_option_count][2][i], label="", linewidth=1, color = "mediumseagreen")
    plt.plot(xdata, ydata_avg_hole_sd_fit,  linewidth=3, color="purple", label="Hole")
    plt.plot(xdata, ydata_avg_elec_sd_fit,  linewidth=3, color="green",  label="Electron")
    plt.tight_layout()

    if dynamics_option == 0:
        plt.savefig("Si75H64_elec_hole_sd_Hot_Energy_Decay_"+decoherence_options[decoherence_option_count])
    elif dynamics_option == 1:
        plt.savefig("Si75H64_elec_hole_sd_Recovery_"+decoherence_options[decoherence_option_count])

    decoherence_option_count += 1


#END
