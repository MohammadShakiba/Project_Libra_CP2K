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
import libra_py.workflows.nbra.lz as lz
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import numpy as np
import matplotlib.pyplot as plt
import os
import multiprocessing as mp




####################
# 1. Get the Hvibs from step2 
print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/_res/")
params["Hvib_re_prefix"] = "Hvib_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_ci_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 4000
params["nstates"]        = 31 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mb = step4.get_Hvib2(params)
hvib_mb[0][-1].show_matrix()
print ("Length of hvib_mb is: ", len(hvib_mb[0]))
#sys.exit(0)

print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/_res/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 4000
params["nstates"]        = 63 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_mixed_sd = step4.get_Hvib2(params)
hvib_mixed_sd[0][-1].show_matrix()
print ("Length of hvib_mixed_sd is: ", len(hvib_mixed_sd[0]))
#sys.exit(0)

print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/res")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 4000
params["nstates"]        = 14 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_elec_sd = step4.get_Hvib2(params)
hvib_elec_sd[0][-1].show_matrix()
print ("Length of hvib_elec_sd is: ", len(hvib_elec_sd[0]))
#sys.exit(0)

print ("\nGathering data from MD ")
absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/../../step3/res/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 4000
params["nstates"]        = 25 # total number of electronic states
params["init_times"]     = [0]
params["active_space"]   = list(range(params["nstates"])) # indexing is from 0!
# Include HOMO and up to the last electronic state
hvib_hole_sd = step4.get_Hvib2(params)
hvib_hole_sd[0][-1].show_matrix()
print ("Length of hvib_hole_sd is: ", len(hvib_hole_sd[0]))
#sys.exit(0)


#####################
# 2. Divide up into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories
subtraj_time_info = [
                      [0, 0], [0, 250], [0, 500], [0, 750], [0, 1000], [0, 1250], [0, 1500], [0,1750], [0, 2000],
                    ]

nsubtrajs   = len(subtraj_time_info)
subtraj_len = 2000 #params["nfiles"] # steps
params["dt"] = 1.0 * units.fs2au

hvib_mb_subtrajs = []
hvib_mixed_sd_subtrajs = []
hvib_elec_sd_subtrajs = []
hvib_hole_sd_subtrajs = []

nstates_mb = hvib_mb[0][0].num_of_rows
nstates_mixed_sd = hvib_mixed_sd[0][0].num_of_rows
nstates_elec_sd  = hvib_elec_sd[0][0].num_of_rows
nstates_hole_sd  = hvib_hole_sd[0][0].num_of_rows

for i in range( nsubtrajs ):
    hvib_mb_subtrajs.append( [] )
    hvib_mixed_sd_subtrajs.append( [] )
    hvib_elec_sd_subtrajs.append( [] )
    hvib_hole_sd_subtrajs.append( [] )


def compute_state_energies_vs_time( hvib ):     
    """
    Computes the states energies vs time for a given hvib.
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
    returns: a list of energies vs time for all states in hvib
    """
    nsteps   = len(hvib) 
    nstates  = hvib[0].num_of_rows
    energies = []
    for state in range( nstates ):
        energies.append( [] )
        for step in range( nsteps ):
            energies[ state ].append( hvib[ step ].get( state, state ).real - hvib[ step ].get( 0, 0 ).real )
    return np.array( energies )


def run_nbra_namd_wrapper( hvib, istate, outfile_name, params ):
    """
    A small wrapper function to run the nbra namd dynamics. 
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
        istate ( int ): Index of the initial state. This index is from 0
        outfile ( string ): The name of the output file
        params ( dict ): A dictionary of important dynamical parameters
    returns: A list of energy values that is the electronic energy vs time. Returned energy is in eV
    """
  
    params["outfile"]  = outfile_name
    params["nstates"]  = hvib[0].num_of_rows
    params["istate"]   = istate
    res_nbra_namd = step4.run( [ hvib ], params )
    nbra_namd     = data_conv.MATRIX2nparray( res_nbra_namd )
    # For SH energies
    energy_nbra_namd = ( nbra_namd[:,3*params["nstates"]+2] - nbra_namd[:,1] )*units.au2ev
    # For SE energies
    #energy_nbra_namd = ( bllz_namd[:,3*params["nstates"]+1] - bllz_namd[:,1] )*units.au2ev
    return energy_nbra_namd



def run_bllz_wrapper( hvib, istate, outfile_name, params ):
    """
    A small wrapper function to run the bllz dynamics. 
        hvib ( list of list of matrices ): The vibronic hamiltonian for all timesteps
        istate ( int ): Index of the initial state. This index is from 0
        outfile ( string ): The name of the output file
        params ( dict ): A dictionary of important dynamical parameters
    returns: A list of energy values that is the electronic energy vs time. Returned energy is in eV
    """

    params["outfile"]  = outfile_name
    params["nstates"]  = hvib[0].num_of_rows
    params["istate"]   = istate
    res_bllz, P = lz.run( [hvib], params )
    bllz_namd   = data_conv.MATRIX2nparray( res_bllz )
    # Recall that for bllz, we use the schrodiger not the surface hopping energy
    #energy_bllz = ( bllz_namd[:,3*params["nstates"]+1] - bllz_namd[:,1] )*units.au2ev
    # for bllz surface hopping energy
    energy_bllz = ( bllz_namd[:,3*params["nstates"]+2] - bllz_namd[:,1] )*units.au2ev
    return energy_bllz



#####################
# 3. Define a function to run the namd
def myfunc( subtraj ):

    nuclear_traj = subtraj_time_info[subtraj][0]
    start_time   = subtraj_time_info[subtraj][1]

    print("\nsubtraj", subtraj)
    print("nuclear trajectory", nuclear_traj)
    print("start time =", start_time)
    print("end time   =", start_time + subtraj_len)

    md_time = [ i for i in range( subtraj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    # Obtain the subtrajectory
    for i in range(start_time, start_time + subtraj_len):
        hvib_mb_subtrajs[ subtraj ].append( hvib_mb[0][i] )
        hvib_mixed_sd_subtrajs[ subtraj ].append( hvib_mixed_sd[0][i] )
        hvib_elec_sd_subtrajs[ subtraj ].append( hvib_elec_sd[0][i] )
        hvib_hole_sd_subtrajs[ subtraj ].append( hvib_hole_sd[0][i] )

    params["nsteps"] = subtraj_len
    hvib_mb_subtraj       = compute_state_energies_vs_time( hvib_mb_subtrajs[ subtraj ] )
    hvib_mixed_sd_subtraj = compute_state_energies_vs_time( hvib_mixed_sd_subtrajs[ subtraj ] )
    hvib_elec_sd_subtraj  = compute_state_energies_vs_time( hvib_elec_sd_subtrajs[ subtraj ] )
    hvib_hole_sd_subtraj  = compute_state_energies_vs_time( hvib_hole_sd_subtrajs[ subtraj ] )

    #========================================================
    # Setting up now for nbra namd calculations
    params["T"]          = 300.0
    params["ntraj"]      = 250
    params["sh_method"]  = 1
    params["Boltz_opt"]  = 1
    params["nsteps"]     = subtraj_len
    params["init_times"] = [0]

    # Many-body, Single-Particle
    istates = [25,26,27,28]
    for istate in istates:
 
        ### FSSH
        params["decoherence_method"] = 0
        outfile = "_out_cd33se33_mb_nbra_namd_fssh_"+str(subtraj)+"_"+str(istate)+".txt"
        mb_nbra_fssh = run_nbra_namd_wrapper( hvib_mb_subtrajs[ subtraj ], istate, outfile, params )
        outfile = "_out_cd33se33_mixed_sd_nbra_namd_fssh_"+str(subtraj)+"_"+str(istate)+".txt"
        mixed_sd_nbra_fssh = run_nbra_namd_wrapper( hvib_mixed_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### ID-A
        params["decoherence_method"] = 1
        outfile = "_out_cd33se33_mb_nbra_namd_ida_"+str(subtraj)+"_"+str(istate)+".txt"
        mb_nbra_ida = run_nbra_namd_wrapper( hvib_mb_subtrajs[ subtraj ], istate, outfile, params )
        outfile = "_out_cd33se33_mixed_sd_nbra_namd_ida_"+str(subtraj)+"_"+str(istate)+".txt"
        mixed_sd_nbra_ida = run_nbra_namd_wrapper( hvib_mixed_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### mSDM
        params["decoherence_method"] = 2
        outfile = "_out_cd33se33_mb_nbra_namd_msdm_"+str(subtraj)+"_"+str(istate)+".txt"
        mb_nbra_msdm = run_nbra_namd_wrapper( hvib_mb_subtrajs[ subtraj ], istate, outfile, params )
        outfile = "_out_cd33se33_mixed_sd_nbra_namd_msdm_"+str(subtraj)+"_"+str(istate)+".txt"
        mixed_sd_nbra_msdm = run_nbra_namd_wrapper( hvib_mixed_sd_subtrajs[ subtraj ], istate, outfile, params )

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ Many-body traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',   fontsize=10)
        plt.ylabel('Energy, eV', fontsize=10)
        plt.ylim(0,3.0)
        plt.yticks([0,1,2,3])
        for mb_index in range( nstates_mb ):
            plt.plot(md_time, hvib_mb_subtraj[mb_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, mb_nbra_fssh, linewidth=2.5, color="black",  label="FSSH")
        plt.plot(md_time, mb_nbra_ida,  linewidth=2.5, color="blue",   label="ID-A")
        plt.plot(md_time, mb_nbra_msdm, linewidth=2.5, color="maroon", label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("Cd33Se33_mb_"+str(subtraj)+"_"+str(istate)+".png", dpi=300)

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ Single-Particle traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',   fontsize=10)
        plt.ylabel('Energy, eV', fontsize=10)
        plt.ylim(0,3.0)
        plt.yticks([0,1,2,3])
        for sp_index in range( nstates_mixed_sd ):
            plt.plot(md_time, hvib_mixed_sd_subtraj[sp_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, mixed_sd_nbra_fssh, linewidth=2.5, color="black",  label="FSSH")
        plt.plot(md_time, mixed_sd_nbra_ida,  linewidth=2.5, color="blue",   label="ID-A")
        plt.plot(md_time, mixed_sd_nbra_msdm, linewidth=2.5, color="maroon", label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("Cd33Se33_mixed_sd_"+str(subtraj)+"_"+str(istate)+".png", dpi=300)

    # Electrons
    istates = [6,7,8,9]
    for istate in istates:

        ### FSSH
        params["decoherence_method"] = 0
        outfile = "_out_cd33se33_elec_sd_nbra_namd_fssh_"+str(subtraj)+"_"+str(istate)+".txt"
        elec_sd_nbra_fssh = run_nbra_namd_wrapper( hvib_elec_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### ID-A
        params["decoherence_method"] = 1
        outfile = "_out_cd33se33_elec_sd_nbra_namd_ida_"+str(subtraj)+"_"+str(istate)+".txt"
        elec_sd_nbra_ida = run_nbra_namd_wrapper( hvib_elec_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### mSDM
        params["decoherence_method"] = 2
        outfile = "_out_cd33se33_elec_sd_nbra_namd_msdm_"+str(subtraj)+"_"+str(istate)+".txt"
        elec_sd_nbra_msdm = run_nbra_namd_wrapper( hvib_elec_sd_subtrajs[ subtraj ], istate, outfile, params )

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ e$^{-}$ traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',   fontsize=10)
        plt.ylabel('Energy, eV', fontsize=10)
        plt.ylim(0,3.0)
        plt.yticks([0,1,2,3])
        for sp_index in range( nstates_elec_sd ):
            plt.plot(md_time, hvib_elec_sd_subtraj[sp_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, elec_sd_nbra_fssh, linewidth=2.5, color="black",  label="FSSH")
        plt.plot(md_time, elec_sd_nbra_ida,  linewidth=2.5, color="blue",   label="ID-A")
        plt.plot(md_time, elec_sd_nbra_msdm, linewidth=2.5, color="maroon", label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("Cd33Se33_elec_sd_"+str(subtraj)+"_"+str(istate)+".png", dpi=300)

    # Hole
    istates = [14,15,16,17,18]
    for istate in istates:

        ### FSSH
        params["decoherence_method"] = 0
        outfile = "_out_cd33se33_hole_sd_nbra_namd_fssh_"+str(subtraj)+"_"+str(istate)+".txt"
        hole_sd_nbra_fssh = run_nbra_namd_wrapper( hvib_hole_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### ID-A
        params["decoherence_method"] = 1
        outfile = "_out_cd33se33_hole_sd_nbra_namd_ida_"+str(subtraj)+"_"+str(istate)+".txt"
        hole_sd_nbra_ida = run_nbra_namd_wrapper( hvib_hole_sd_subtrajs[ subtraj ], istate, outfile, params )

        ### mSDM
        params["decoherence_method"] = 2
        outfile = "_out_cd33se33_hole_sd_nbra_namd_msdm_"+str(subtraj)+"_"+str(istate)+".txt"
        hole_sd_nbra_msdm = run_nbra_namd_wrapper( hvib_hole_sd_subtrajs[ subtraj ], istate, outfile, params )

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ h$^{+}$ traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',   fontsize=10)
        plt.ylabel('Energy, eV', fontsize=10)
        plt.ylim(0,3.0)
        plt.yticks([0,1,2,3])
        for sp_index in range( nstates_hole_sd ):
            plt.plot(md_time, hvib_hole_sd_subtraj[sp_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, hole_sd_nbra_fssh, linewidth=2.5, color="black",  label="FSSH")
        plt.plot(md_time, hole_sd_nbra_ida,  linewidth=2.5, color="blue",   label="ID-A")
        plt.plot(md_time, hole_sd_nbra_msdm, linewidth=2.5, color="maroon", label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("Cd33Se33_hole_sd_"+str(subtraj)+"_"+str(istate)+".png", dpi=300)

    #start_time += subtraj_increment




pool = mp.Pool(nsubtrajs)
pool.map( myfunc, list(range(nsubtrajs)) )
pool.close()
pool.join()




