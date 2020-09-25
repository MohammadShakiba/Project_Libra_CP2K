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

colors = {}
colors.update({"1": '#000000'})  # Black 
colors.update({"2": '#000099'})  # Blue  
colors.update({"3": '#006400'})  # Green 
colors.update({"4": '#990000'})  # Red   
colors.update({"5": '#8B008B'})  # Purple
colors.update({"6": '#FF8C00'})  # Orange
colors.update({"9": '#4d4d4d'})  # Gray  
color_index = ["1","2","3","4","5","6","9"]




# 1) Define the fitting function(s).
def func3(t, tau, beta, E0, B):
    return E0*np.exp( -(t/tau)**beta ) + B



# 2) Set paramters for reading / sorting the data
dt = 1.0 #fs

basis_options = ["sd","mb"] 
decoherence_options = ["fssh","ida", "msdm"]#, "dish"]
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
istates = [ [33,34,35], [33,34,35] ]

r_squared_tol = 0.8

taus = []
dyns = []

decoherence_option_count = 0
# 1. For each decoherence option, we will perofrm dynamics using both the SD and MB bases
for decoherence_option in decoherence_options:

    print("Number of initial conditions before R2 filtering", nsubtrajs*len(istates[0]))

    taus.append( [] )
    dyns.append( [] )

    # 2. For each basis option choice (SD and MB) we have a number of initial conditions, we will obtain a tau for each
    basis_option_count = 0
    for basis_option in basis_options:
        
        taus[decoherence_option_count].append( [] )
        dyns[decoherence_option_count].append( [] )     

        # 3. Now, the initial conditions are for 4 nuclear sub-trajectories, each having 3 initial "istate" conditions
        for subtraj in range(nsubtrajs):

            # 3.1. Again, each istate is also an initial condition, here. So for each initial condition, we collect a tau
            for istate in istates[basis_option_count]:

                xdata     = []
                ydata     = []
                ydata_fit = []

                filename = "../run_namd/_out_"+basis_option+"_"+decoherence_option+"_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
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
                        # sp dynamics
                        if basis_option_count == 0:
                            y = float( namd_data_line[ 194 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.0001

                        # mb dynamics
                        else:
                            y = float( namd_data_line[ 95 ] ) - float( namd_data_line[4] )
                            if y < 0:
                                y = 0.0001
                        y = y * units.au2ev
                        #print(y)

                    # Population Recovery Dynamics
                    elif dynamics_option == 1:
                        if basis_option_count == 0 or basis_option_count == 1:
                            y = 1.0 - ( float( namd_data_line[ 6 ] ) )

                    # Ground state recovery
                    elif dynamics_option == 2:
                        if basis_option_count == 0 or basis_option_count == 1:
                            y = 1.0 - ( float( namd_data_line[ 3 ] ) ) 

                    xdata.append(i*dt)
                    ydata.append(y)

                dyns[decoherence_option_count][basis_option_count].append( np.array( ydata ) )

                # At this point, we have obtained the x and y data for this particular combination of dynamics, basis, and decoherence option.
                # Now, to fit the data 
                ydata_fit = []

                #if decoherence_options[decoherence_option_count] == "msdm" and basis_option == "sd" and subtraj == 9 and istate == 6:
                #    print("BREAKBREAK")
                #    break 
                #if decoherence_options[decoherence_option_count] == "msdm" and basis_option == "mb" and subtraj == 20 and istate == 6:
                #    print("BREAKBREAK")
                #    break

                print ("\nFitting with the stretched-compressed exponential fitting function" )
                popt, pcov = curve_fit(func3, xdata, ydata, bounds=([0.0, 0.0, ydata[0]-0.001, 0.0], [np.inf, np.inf, ydata[0]+0.001, 0.001]))
                residuals  = ydata - func3(xdata, *popt)
                ss_res     = np.sum(residuals**2)
                ss_tot     = np.sum((ydata - np.mean(ydata))**2)
                r_squared  = 1.0 - (ss_res / ss_tot)
                #print ("r_squared = ", r_squared)
                tau, beta, E0, B = popt
                #print ("Printing optimal fitting paramters:")
                #print ("B    = ", B)
                #print ("E0   = ", E0, "eV")
                #print ("beta = ", beta)
                #print ("tau  = ", tau, "fs")
                for i in range(len(xdata)):
                    ydata_fit.append(func3(xdata[i], *popt))          
 
                plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
                plt.subplot(1,1,1)
                plt.title("Si$_{75}$H$_{64}$ traj"+str(subtraj)+" "+basis_option+" "+decoherence_option+" istate = "+str(istate), fontsize=8) 
                plt.xlabel('Time, fs',   fontsize=8)
                plt.ylabel('Population', fontsize=8)
                plt.yticks([0.0,0.2,0.4,0.6,0.8,1.0])
                #plt.ylabel('Hot Energy, eV',  fontsize=8)
                #plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7])
                plt.plot(xdata, ydata,     linewidth=1.5, color="black", label="data")
                plt.plot(xdata, ydata_fit, linewidth=1,   color="green", label="fit") 
                plt.legend(fontsize=8)
                plt.tight_layout()
                plt.savefig(outname)
                #sys.exit(0)

                # Append a tau for each initial condition, here, initial conditions are istates ..
                if r_squared > r_squared_tol:
                    taus[decoherence_option_count][basis_option_count].append(tau)
                    print ("r_squared = ", r_squared)
                    print ("Printing optimal fitting paramters:")
                    print ("B    = ", B)
                    print ("E0   = ", E0, "eV")
                    print ("beta = ", beta)
                    print ("tau  = ", tau, "fs")


        niconds = len(taus[decoherence_option_count][basis_option_count])
        print("final ", decoherence_options_names[decoherence_option_count], basis_options[basis_option_count], "niconds", niconds)

        taus_avg, taus_std = data_stat.scalar_stat( taus[decoherence_option_count][basis_option_count] )       
        error_95 = 1.96 * taus_std / math.sqrt(niconds)

        print ("\nFinished computing data for ", decoherence_options[decoherence_option_count], [basis_option_count])
        print ("Average decay lifetime over every NAD run: ",             taus_avg, " fs")
        print ("Std_dev decay lifetime over every NAD run: ",             taus_std, " fs")
        print ("Average population decay lifetime with error (95% ci): ", taus_avg, "+-", error_95, " fs\n")
        print ("final", decoherence_options[decoherence_option_count], basis_options[basis_option_count],  taus_avg, "+-", error_95, " fs\n") 

        basis_option_count += 1

    # 4. We have compute the dynamics for all bases (SD and/or MB) for a given decoherence scheme. Now, let's make our plots
    ydata_avg_sd = sum(dyns[decoherence_option_count][0]) / len(dyns[decoherence_option_count][0])
    ydata_avg_sd_fit = []

    ydata_avg_mb = sum(dyns[decoherence_option_count][1]) / len(dyns[decoherence_option_count][1])
    ydata_avg_mb_fit = []

    print ("\nFitting with the stretched-compressed exponential fitting function" )
    popt, pcov = curve_fit(func3, xdata, ydata_avg_mb, bounds=([0.0, 0.0, ydata_avg_mb[0]-0.001, 0.0], [np.inf, np.inf, ydata_avg_mb[0]+0.001, 0.01]))
    residuals  = ydata_avg_mb - func3(xdata, *popt)
    ss_res     = np.sum(residuals**2)
    ss_tot     = np.sum((ydata_avg_mb - np.mean(ydata_avg_mb))**2)
    r_squared  = 1.0 - (ss_res / ss_tot)
    print ("r_squared = ", r_squared)
    tau, beta, E0, B = popt
    print ("Printing optimal fitting paramters:")
    print ("B    = ", B)
    print ("E0   = ", E0, "eV")
    print ("beta = ", beta)
    print ("tau  = ", tau, "fs")
    for i in range(len(xdata)):
        ydata_avg_mb_fit.append(func3(xdata[i], *popt))
    print("avg mb tau = ", tau)

    print ("\nFitting with the stretched-compressed exponential fitting function" )
    popt, pcov = curve_fit(func3, xdata, ydata_avg_sd, bounds=([0.0, 0.0, ydata_avg_sd[0]-0.001, 0.0], [np.inf, np.inf, ydata_avg_sd[0]+0.001, 0.01]))
    residuals  = ydata_avg_sd - func3(xdata, *popt)
    ss_res     = np.sum(residuals**2)
    ss_tot     = np.sum((ydata_avg_sd - np.mean(ydata_avg_sd))**2)
    r_squared  = 1.0 - (ss_res / ss_tot)
    print ("r_squared = ", r_squared)
    tau, beta, E0, B = popt
    print ("Printing optimal fitting paramters:")
    print ("B    = ", B)
    print ("E0   = ", E0, "eV")
    print ("beta = ", beta)
    print ("tau  = ", tau, "fs")
    for i in range(len(xdata)):
        ydata_avg_sd_fit.append(func3(xdata[i], *popt))
    print("avg sd tau = ", tau)

    plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
    plt.subplot(1,1,1)
    plt.title("Si$_{75}$H$_{64}$"+decoherence_options_names[decoherence_option_count], fontsize=12)
    plt.xlabel('Time, fs',        fontsize=12)
    plt.xticks(fontsize=12) 

    if dynamics_option == 0:
        plt.ylabel('Hot Energy, eV', fontsize=12)
        plt.yticks([0.0,0.2,0.4,0.6,0.8],fontsize=12)
        plt.ylim(0,0.8)
    elif dynamics_option == 1:
        plt.ylabel('Population',     fontsize=12)
        #plt.yticks([0.0,0.2,0.4,0.6,0.8,1.0],fontsize=12)
        plt.yticks([0.0,0.3,0.6,1.0],fontsize=12)
        plt.ylim(0,1.0)

    for i in range( niconds ):
        plt.plot( xdata, dyns[decoherence_option_count][0][i], label="", linewidth=1, color = "pink")
        plt.plot( xdata, dyns[decoherence_option_count][1][i], label="", linewidth=1, color = "skyblue")
    plt.plot(xdata, ydata_avg_sd,     linewidth=2.5, color="black", label="")
    plt.plot(xdata, ydata_avg_mb,     linewidth=2.5, color="black", label="")
    plt.plot(xdata, ydata_avg_sd_fit, linewidth=1,   color="red",   label="Single-particle")
    plt.plot(xdata, ydata_avg_mb_fit, linewidth=1,   color="blue",  label="Many-body")

    #if decoherence_option_count == 0 and dynamics_option == 1:
    #    plt.legend(fontsize=9, ncol=2, loc="upper right")
    #else:
    #    pass
    plt.tight_layout()

    if dynamics_option == 0:
        plt.savefig("Si75H64_Hot_Energy_Decay_"+decoherence_options[decoherence_option_count])
    elif dynamics_option == 1:
        plt.savefig("Si75H64_Recovery_"+decoherence_options[decoherence_option_count])

    decoherence_option_count += 1

