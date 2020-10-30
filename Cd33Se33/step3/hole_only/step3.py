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
data_dim = 82 # rows in E_ks
active_space = range(data_dim)
start_time   = 0   
finish_time  = 2#4001   
dt = 1.0*units.fs2au
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
St_ks[0][-2].show_matrix()
#sys.exit(0)




####################
# 2.0. Generate the homo-n -> lumo excitations. We know before hand that homo-23 -> lumo is to be the highest considered excitation. Spin-polarization is not considered in this study, so excitation of only 1 spin type (alpha herein) is chosen. These "excitations" are in the format of how CP2K would output them. We do this because the routine used below can transform this format of "excitations" into a format expected by Libra
unique_SDs = []
unique_SDs.append(['28, 28'])
for i in range(28,4,-1):
    unique_SDs.append(['%d, 29'%i])
print(unique_SDs)
#sys.exit(0)





####################
#2.1. Here, we transform the unique SD excitatons as given by CP2K into the format Libra expects. We considered only transitions within the spin-restricted reference, so our excitations are of 1 spin type. Here we choose the alpha electrons only. 
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
# 4.0 Compute the mid-point energies and time-overlaps according to the energy ordered SDs.
St_sd = []
E_sd  = []
res_dir_sd = "res_sd"
os.system("rm -r "+res_dir_sd )
os.system("mkdir "+res_dir_sd )
for step in range(finish_time-start_time-1):
    E_sd.append ( 0.5 * (E_sd_sorted[step] + E_sd_sorted[step+1]) )
    St_sd.append( mapping.ovlp_mat_arb( basis_sorted[step], basis_sorted[step+1], St_ks[0][step], use_minimal=False ) )




####################
# 5.0. Apply state reordering and/or phase corrections to the St_sd data
params["do_state_reordering"]    = 2
params["state_reordering_alpha"] = 0.0
step3.apply_state_reordering_general(St_sd, E_sd, params)
step3.apply_phase_correction_general( St_sd )

# 4. Form Hvib_sd
Hvib_sd = []
for step in range(finish_time-start_time-1):
    hvib_sd = E_sd[step] - 0.5j/dt * (St_sd[step] - St_sd[step].H())
    Hvib_sd.append( hvib_sd )
    Hvib_sd[step].real().show_matrix("%s/Hvib_sd_%d_re" % (res_dir_sd, step))
    Hvib_sd[step].imag().show_matrix("%s/Hvib_sd_%d_im" % (res_dir_sd, step))
    St_sd[step].real().show_matrix("%s/St_%d_re" % (res_dir_sd, step))
#sys.exit(0)

