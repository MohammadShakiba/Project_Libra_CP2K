import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize
from libra_py import data_conv
from libra_py import data_read
from libra_py import data_stat
import libra_py.workflows.nbra.mapping as mapping
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.step3 as step3
import numpy as np
import matplotlib.pyplot as plt




###############
# 1. Read the files that have the energies and time-overlap matricies in the Kohn-Sham basis, E_ks and St_ks.
path = os.getcwd()
res_dir = path+"/../../step2/res/"
data_dim = 78 # rows in E_ks
active_space = range(data_dim)
start_time   = 2000     # initial step
finish_time  = 6001  # final step + 1  
dt = 0.5*units.fs2au
params = { "data_set_paths" : [res_dir], "data_dim":data_dim, "active_space":active_space, "isnap":start_time,  "fsnap":finish_time }
# Fetching E_ks
params.update({ "data_re_prefix" : "E_ks_",  "data_re_suffix" : "_re", "data_im_prefix" : "E_ks_",  "data_im_suffix" : "_im"  } )
E_ks = data_read.get_data_sets(params)
print ("\n")
E_ks[0][-1].show_matrix()
#sys.exit(0)
# Fetching St_ks
params.update({ "data_re_prefix" : "St_ks_", "data_re_suffix" : "_re", "data_im_prefix" : "St_ks_", "data_im_suffix" : "_im"  } )
St_ks = data_read.get_data_sets(params)
print ("\n")
St_ks[0][-1].show_matrix()
#sys.exit(0)





####################
# 2.0. shown below are excitations. These are obtained from running "get_unique.py", which goes into the slurm output file (that is, output by Libra) that was made from processing the output files of the TD-DFT cp2k calculations and extracts the unique set of SDs over all timesteps. These SDs are given in terms of their 1-electron KS excitation from a homo index w.r.t 1. We then manually added the ground state as ['31, 31'], where 31 is the homo index from 1. This show some excitaiton from 31 to itself, which results in no change, and is thus the HOMO. See section 2.1.

unique_SDs = [ ['31, 31'], ['31, 32'], ['30, 32'], ['29, 32'], ['31, 33'], ['28, 32'], ['30, 33'], ['29, 33'], ['27, 32'], ['28, 33'], ['26, 32'], ['27, 33'], ['25, 32'], ['31, 34'], ['26, 33'], ['24, 32'], ['30, 34'], ['25, 33'], ['31, 35'], ['23, 32'], ['29, 34'], ['24, 33'], ['28, 34'], ['30, 35'], ['29, 35'], ['23, 33'], ['22, 32'], ['27, 34'], ['21, 32'], ['28, 35'], ['26, 34'], ['31, 36'], ['20, 32'], ['22, 33'], ['27, 35'], ['21, 33'], ['31, 37'], ['26, 35'], ['25, 34'], ['29, 36'], ['19, 32'], ['30, 36'], ['24, 34'], ['20, 33'], ['25, 35'], ['30, 37'], ['23, 34'], ['19, 33'], ['18, 32'], ['29, 37'], ['28, 36'], ['24, 35'], ['17, 32'], ['31, 38'], ['27, 36'], ['30, 38'], ['31, 39'], ['31, 40'], ['30, 39'], ['31, 41'], ['28, 37'], ['22, 34'], ['29, 38'], ['28, 38'], ['30, 40'], ['29, 39'], ['16, 32'], ['23, 35'], ['26, 36'], ['31, 42'], ['21, 34'], ['25, 36'], ['15, 32'], ['27, 37'], ['14, 32'], ['13, 32'], ['12, 32'], ['11, 32']]

# 2.1. Here, we transform the unique SD excitatons as given by CP2K into the format Libra expects. We considered only transitions within the spin-restricted reference, so our excitations are of 1 spin type. Here we choose the alpha electrons only. 
max_n = 0
for i in range(len(unique_SDs)):
    max_n = max([int(unique_SDs[i][0].split(',')[0]), max_n])
basis = []
for i in range(len(unique_SDs)):
    tmp_basis = [x for x in range(1,max_n+1)]
    tmp_basis[int(unique_SDs[i][0].split(',')[0])-1] = int(unique_SDs[i][0].split(',')[1])
    basis.append(tmp_basis)
print (basis)




####################
# 3.0. Now, for each timestep, we compute the energies of the SDs and sort them by energy. This prevents the state energies from crossing in this basis. However, note that the St_ks files (time-overlaps in the KS basis) have no knowledge of the reordering of the SDs by energies
nSDs = len(basis)
E_sd_sorted  = []
basis_sorted = []
SD_energy_corr = [0.0]*nSDs
for step in range(finish_time-start_time):

    E = mapping.energy_mat_arb( basis, E_ks[0][step], SD_energy_corr )
    nstates = E.num_of_cols
    e = np.zeros( nstates )
    for state in range(nstates):
        e[state] =  E.get(state,state).real
    reindex = np.argsort(e)
    E_sd_sorted.append( CMATRIX(nstates,nstates) )
    basis_sorted.append( [] )
    for i in range(len(reindex)):
        E_sd_sorted[step].set( i, i, E.get(  int(reindex[i]), int(reindex[i]) ) )
        basis_sorted[step].append( basis[ int(reindex[i]) ] )



####################
# 4.0 Compute the mid-point energies and time-overlaps according to the energy ordered SDs. Note that by computing wavefunciton overlaps of the energy ordered SDs using the properties of the KS orbtials results in the time-overlaps in the SD basis taking their largest magnitudes in off-diagonal locations. This is because the ordering of the wavefuncitons in the KS basis is not consistent with the ordering of the SD basis based on energy ordering. Only in special cases does the energy ordering of the SD basis match the ordering of the KS basis (such as for homo->lumo+N or homo-N->lumo only excitations). Therefore, we will need to apply a state reordering procedure to the St_sd matricies computed here.
St_sd = []
E_sd  = []
res_dir_sd = "res_sd"
os.system("rm -r "+res_dir_sd )
os.system("mkdir "+res_dir_sd )
for step in range(finish_time-start_time-1):
    E_sd.append ( 0.5 * (E_sd_sorted[step] + E_sd_sorted[step+1]) )
    St_sd.append( mapping.ovlp_mat_arb( basis_sorted[step], basis_sorted[step+1], St_ks[0][step], use_minimal=False ) )




####################
# 5.0 Now, apply the state reordering procedured mentioned above, followed by a correction to the phases.
params["do_state_reordering"]    = 2
params["state_reordering_alpha"] = 0.0
step3.apply_state_reordering_general(St_sd, E_sd, params)
step3.apply_phase_correction_general( St_sd )




####################
# 6.0 Form Hvib_sd
Hvib_sd = []
for step in range(finish_time-start_time-1):
    hvib_sd = E_sd[step] + 0.5j/dt * (St_sd[step] - St_sd[step].H())
    Hvib_sd.append( hvib_sd )
    Hvib_sd[step].real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir_sd, step))
    Hvib_sd[step].imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir_sd, step))
    St_sd[step].real().show_matrix("%s/St_%d_re" % (res_dir_sd, step))

