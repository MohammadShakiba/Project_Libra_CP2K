import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_stat
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import numpy as np
import matplotlib.pyplot as plt
import os

colors = {}
colors.update({"1": '#000000'})  # Black 
colors.update({"2": '#000099'})  # Blue  
colors.update({"3": '#006400'})  # Green 
colors.update({"4": '#990000'})  # Red   
colors.update({"5": '#8B008B'})  # Purple
colors.update({"6": '#FF8C00'})  # Orange
colors.update({"9": '#4d4d4d'})  # Gray  
color_index = ["1","2","3","4","5","6","9"]


os.system("rm _*")
print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/mixed_electron_hole/")
params["Hvib_re_prefix"] = "Hvib_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_ci_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 4000
params["nstates"]        = 31 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list( range(params["nstates"]) ) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mb = step4.get_Hvib2(params)
hvib_mb[0][-1].show_matrix()
print ("Length of hvib_mb is: ", len(hvib_mb[0]))
#sys.exit(0)

params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/mixed_electron_hole/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nstates"]        = 63 # total number of electronic states
params["nfiles"]         = 4000
params["init_times"]     = [0]
params["active_space"]   = list( range(params["nstates"]) ) # indexing is from 0!
hvib_sd = step4.get_Hvib2(params)
hvib_sd[0][-1].show_matrix()
print ("Length of hvib_sd is: ", len(hvib_sd[0]))
#sys.exit(0)





# 2. Divide up into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories
nsubtrajs   = 1
subtraj_len = params["nfiles"] # steps, so 500 fs with 0.5 dt
params["dt"] = 1.0 * units.fs2au

hvib_mb_subtrajs = []
hvib_sd_subtrajs = []

for subtraj in range(nsubtrajs):

    hvib_mb_subtrajs.append( [] )
    hvib_sd_subtrajs.append( [] )

    start_time = subtraj*subtraj_len
    print("\nstart time = ", start_time)
    print("end time   = ", start_time + subtraj_len)

    md_time = [ i for i in range( subtraj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    # Obtain the subtrajectory
    for i in range(start_time, start_time + subtraj_len):
        hvib_mb_subtrajs[ subtraj ].append( hvib_mb[0][i] )
        hvib_sd_subtrajs[ subtraj ].append( hvib_sd[0][i] )


    # Obtain all nacs for all pairs states 
    nac_mb, nac_sd = [], []
    states_mb = range( hvib_mb[0][0].num_of_rows )
    states_sd = range( hvib_sd[0][0].num_of_rows )

    nac_cutoff = 0 #meV
    mvalues    = [0]

    for mvalue in range(len(mvalues)):

        nac_mb.append( [] )
        nac_sd.append( [] )

        for t in range(start_time, start_time + subtraj_len):

            # MB
            for i in states_mb:
                for j in states_mb:    
                    if j != i:

                        nac_mb_ij = abs( hvib_mb[0][t].get(i,j).imag ) * 1000.0 * units.au2ev
                        if nac_mb_ij > nac_cutoff:
                            x_mb = MATRIX(1,1)
                            x_mb.set(0, 0, nac_mb_ij )
                            nac_mb[mvalue].append( x_mb )

            # SD
            for i in states_sd:
                for j in states_sd:
                    if j != i:
                        nac_sd_ij = abs( hvib_sd[0][t].get(i,j).imag ) * 1000.0 * units.au2ev
                        if nac_sd_ij > nac_cutoff:
                            x_sd = MATRIX(1,1)
                            x_sd.set(0, 0, nac_sd_ij )
                            nac_sd[mvalue].append( x_sd )


        #==========================
        # Sanity check
        print("len(hvib_mb[0]) is :",len(hvib_mb[0]))
        nac_mb[0][0].show_matrix()
        nac_mb[0][-1].show_matrix()
        print( hvib_mb[0][0].get(0,1)    * 1000.0 * units.au2ev )
        print( hvib_mb[0][-1].get(29,30) * 1000.0 * units.au2ev )
        print("len(nac_mb[0]) is :",len(nac_mb[0]))
        print("len(nac_sd[0]) is :",len(nac_sd[0]))
        # End of Sanity check
        #==========================

        # For MB
        bin_supp_mb, dens_mb, cum_mb = data_stat.cmat_distrib( nac_mb[mvalue], 0, 0, 0, 0, 150, 0.5)

        # For SD
        bin_supp_sd, dens_sd, cum_sd = data_stat.cmat_distrib( nac_sd[mvalue], 0, 0, 0, 0, 150, 0.5)



fs = 15
plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
plt.subplot(1,1,1)
plt.title('(CdSe)$_{33}$', fontsize=fs)
plt.xlabel('|NAC|, meV',  fontsize=fs)
plt.ylabel('PD, 1/meV', fontsize=fs)
plt.yticks(fontsize=fs)
plt.xticks(fontsize=fs)
plt.xlim(0,10)
plt.plot( bin_supp_sd, dens_sd, label="SP", linewidth=2.0, color = "red")
plt.plot( bin_supp_mb, dens_mb, label="MB",       linewidth=2.0, color = "blue")
plt.legend(fontsize=12, loc='upper right')
plt.tight_layout()
plt.savefig("CdSe_NAC_Dist_MB_vs_SD_0_10.png", dpi=600)

plt.xlim(10,50)
plt.ylim(0,0.015)
plt.yticks(ticks=[0,0.005,0.01,0.015],labels=['0','5','10','15'],fontsize=fs)
plt.text(10,0.0154,'x10$^{-3}$',fontsize=12)
plt.legend('',frameon=False)
plt.tight_layout()
plt.savefig("CdSe_NAC_Dist_MB_vs_SD_10_50.png", dpi=600)

plt.xlim(50,151)
plt.ylim(0,0.0009)
plt.yticks([0,0.0003,0.0006,0.0009,],labels=['0','3','6','9'],fontsize=fs)
plt.text(49,0.00092,'x10$^{-4}$',fontsize=12)
plt.xticks(list(range(50,151,25)))
plt.tight_layout()
plt.savefig("CdSe_NAC_Dist_MB_vs_SD_50_150.png", dpi=600)


