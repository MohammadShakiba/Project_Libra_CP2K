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

import multiprocessing as mp

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
params["data_set_paths"].append(absolute_path+"/../../../step2/res/")
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
params["data_set_paths"].append(absolute_path+"/../../../step3/res_sd/")
params["Hvib_re_prefix"] = "Hvib_sd_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_sd_"; params["Hvib_im_suffix"] = "_im"
params["nstates"]        = 64 # total number of electronic states
params["nfiles"]         = 4000
params["init_times"]     = [0]
params["active_space"]   = list( range(params["nstates"]) ) # indexing is from 0!
hvib_sd = step4.get_Hvib2(params)
hvib_sd[0][-1].show_matrix()
print ("Length of hvib_sd is: ", len(hvib_sd[0]))
#sys.exit(0)

# 2. Divide up into sub-trajectories. These are to be consdiered our independent nuclear sub-trajectories
nsubtrajs   = 21
subtraj_len = 2000 # steps, so 500 fs with 0.5 dt
start_time = 0
subtraj_increment = 100 # fs
params["dt"] = 1.0 * units.fs2au

hvib_mb_subtrajs = []
hvib_sd_subtrajs = []
for i in range( nsubtrajs ):
    hvib_mb_subtrajs.append( [] )
    hvib_sd_subtrajs.append( [] )

tmp_mb_energies = []
tmp_sd_energies = []

tmp_fssh_mb = []
tmp_ida_mb  = []
tmp_msdm_mb = []

tmp_fssh_sd = []
tmp_ida_sd  = []
tmp_msdm_sd = []


#for subtraj in range(nsubtrajs):

def myfunc( subtraj ):

    
    start_time = subtraj*subtraj_increment

    print("\nsubtraj ", subtraj)
    print("start time = ", start_time)
    print("end time   = ", start_time + subtraj_len)

    md_time = [ i for i in range( subtraj_len ) ]
    md_time = np.array( md_time ) * params["dt"] * units.au2fs

    # Obtain the subtrajectory
    for i in range(start_time, start_time + subtraj_len):
        hvib_mb_subtrajs[ subtraj ].append( hvib_mb[0][i] )
        hvib_sd_subtrajs[ subtraj ].append( hvib_sd[0][i] )

    # Compute mb energy gaps and decoherence times 
    params["nsteps"]     = subtraj_len
    params["init_times"] = [0]
    tau_mb, rates_mb = decoherence_times.decoherence_times( hvib_mb_subtrajs[ subtraj ] )
    dE_mb            = decoherence_times.energy_gaps(       hvib_mb_subtrajs[ subtraj ] )
    avg_deco_mb      = tau_mb * units.au2fs
    print ("Finished gather data for MD time", params["init_times"][0] + params["nsteps"], "steps.")
    print ("Dephasing time between mb states 0 and 1 is:", tau_mb.get(0,1) * units.au2fs, "fs")

    # Compute sd energy gaps and decoherence times 
    params["nsteps"]     = subtraj_len
    params["init_times"] = [0]
    tau_sd, rates_sd = decoherence_times.decoherence_times( hvib_sd_subtrajs[ subtraj ] )
    dE_sd            = decoherence_times.energy_gaps(       hvib_sd_subtrajs[ subtraj ] )
    avg_deco_sd      = tau_sd * units.au2fs
    print ("Finished gather data for MD time", params["init_times"][0] + params["nsteps"], "steps.")
    print ("Dephasing time between sd states 0 and 1 is:", tau_sd.get(0,1) * units.au2fs, "fs")

    # Compute mb energies
    mb_subtraj_energy = []
    params["nstates"] = 31
    for mb_index in range( params["nstates"] ):
        mb_subtraj_energy.append( [] )
        for step in range( params["nsteps"] ):
            mb_subtraj_energy[ mb_index ].append( hvib_mb_subtrajs[ subtraj ][ step ].get( mb_index, mb_index ).real - hvib_mb_subtrajs[ subtraj ][ step ].get( 0, 0 ).real )
    mb_subtraj_energy = np.array( mb_subtraj_energy )
    tmp_mb_energies.append( mb_subtraj_energy )

    # Compute sd energies
    params["nstates"] = 64
    sd_subtraj_energy = []
    for sd_index in range( params["nstates"] ):
        sd_subtraj_energy.append( [] )
        for step in range( params["nsteps"] ):
            sd_subtraj_energy[ sd_index ].append( hvib_sd_subtrajs[ subtraj ][ step ].get( sd_index, sd_index ).real - hvib_sd_subtrajs[ subtraj ][ step ].get( 0, 0 ).real )
    sd_subtraj_energy = np.array( sd_subtraj_energy )
    tmp_sd_energies.append( sd_subtraj_energy )

    ###############################################
    ###############################################
    # Now it is time to run the NAMD
    params["T"]                  = 300.0
    params["ntraj"]              = 250
    params["sh_method"]          = 1
    params["Boltz_opt"]          = 1

    params["nsteps"] = subtraj_len
    params["init_times"] = [0]

    istates = [14,18]

    for istate in istates:

        params["istate"]   = istate # Recall index from 0

        ### FSSH
        params["decoherence_method"] = 0
        start = time.time()
        # mb fssh
        params["outfile"]  = "_out_mb_fssh_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"]  = 31 # total number of mb electronic states
        res_mb_fssh = step4.run( [ hvib_mb_subtrajs[ subtraj ] ], params)
        mb_nbra_namd = data_conv.MATRIX2nparray( res_mb_fssh )
        tmp_fssh_mb.append( mb_nbra_namd )
        mb_hot_energy_fssh = ( mb_nbra_namd[:,95] - mb_nbra_namd[:,1] ) * units.au2ev
        #sys.exit(0)

        #sd fssh
        params["outfile"]  = "_out_sd_fssh_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"]  = 64  # total number of sd electronic states
        res_sd_fssh = step4.run( [ hvib_sd_subtrajs[ subtraj ] ], params)
        sd_nbra_namd  = data_conv.MATRIX2nparray( res_sd_fssh )
        tmp_fssh_sd.append( sd_nbra_namd )
        sd_hot_energy_fssh = ( sd_nbra_namd[:,194] - sd_nbra_namd[:,1] ) * units.au2ev
        end = time.time()
        print("Time to run NBRA NAMD = ", end - start)

        ### IDA
        params["decoherence_method"] = 1
        start = time.time()
        # mb ida
        params["outfile"] = "_out_mb_ida_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"] = 31 # total number of mb electronic states
        res_mb_ida = step4.run( [ hvib_mb_subtrajs[ subtraj ] ], params)
        mb_nbra_namd = data_conv.MATRIX2nparray( res_mb_ida )
        tmp_ida_mb.append( mb_nbra_namd )
        mb_hot_energy_ida = ( mb_nbra_namd[:,95] - mb_nbra_namd[:,1] ) * units.au2ev

        #sd ida
        params["outfile"] = "_out_sd_ida_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"] = 64 # total number of sd electronic states
        res_sd_ida = step4.run( [ hvib_sd_subtrajs[ subtraj ] ], params)
        sd_nbra_namd = data_conv.MATRIX2nparray( res_sd_ida )
        tmp_ida_sd.append( sd_nbra_namd )
        sd_hot_energy_ida = ( sd_nbra_namd[:,194] - sd_nbra_namd[:,1] ) * units.au2ev
        end = time.time()
        print("Time to run NBRA NAMD = ", end - start)

        ### mSDM
        params["decoherence_method"] = 2
        start = time.time()
        # mb msdm
        params["outfile"] = "_out_mb_msdm_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"] = 31 # total number of mb electronic states
        res_mb_msdm = step4.run( [ hvib_mb_subtrajs[ subtraj ] ], params)
        mb_nbra_namd = data_conv.MATRIX2nparray( res_mb_msdm )
        tmp_msdm_mb.append( mb_nbra_namd )
        mb_hot_energy_msdm = ( mb_nbra_namd[:,95] - mb_nbra_namd[:,1] ) * units.au2ev

        #sd msdm
        params["outfile"] = "_out_sd_msdm_subtraj_"+str(subtraj)+"_istate"+str(istate)+".txt"
        params["nstates"] = 64 # total number of sd electronic states
        res_sd_msdm = step4.run( [ hvib_sd_subtrajs[ subtraj ] ], params)
        sd_nbra_namd  = data_conv.MATRIX2nparray( res_sd_msdm )
        tmp_msdm_sd.append( sd_nbra_namd )
        sd_hot_energy_msdm = ( sd_nbra_namd[:,194] - sd_nbra_namd[:,1] ) * units.au2ev
        end = time.time()
        print("Time to run NBRA NAMD = ", end - start)

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ MB traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',      fontsize=10)
        plt.ylabel('Energy, eV',    fontsize=10)
        plt.ylim(0,3.5)
        plt.yticks([0,1,2,3])
        for sd_index in range( 31 ):
            plt.plot(md_time, mb_subtraj_energy[sd_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, mb_hot_energy_fssh, linewidth=2.5, color="red",   label="FSSH")
        plt.plot(md_time, mb_hot_energy_ida,  linewidth=2.5, color="green", label="IDA")
        plt.plot(md_time, mb_hot_energy_msdm, linewidth=2.5, color="blue",  label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("CdSe33_MB_Dyn_"+str(subtraj)+"_"+str(istate)+".png", dpi=300)

        plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)
        plt.subplot(1,1,1)
        plt.title('(CdSe)$_{33}$ SP traj'+str(subtraj)+' istate'+str(istate)+'', fontsize=10)
        plt.xlabel('Time, fs',      fontsize=10)
        plt.ylabel('Energy, eV',    fontsize=10)
        plt.ylim(0,3.5)
        plt.yticks([0,1,2,3])
        for sd_index in range( 64 ):
            plt.plot(md_time, sd_subtraj_energy[sd_index]*units.au2ev, label="", linewidth=1, color = "gray")
        plt.plot(md_time, sd_hot_energy_fssh, linewidth=2.5, color="red",   label="FSSH")
        plt.plot(md_time, sd_hot_energy_ida,  linewidth=2.5, color="green", label="IDA")
        plt.plot(md_time, sd_hot_energy_msdm, linewidth=2.5, color="blue",  label="mSDM")
        plt.tight_layout()
        plt.legend(fontsize=8, ncol=3, loc="lower left")
        plt.savefig("CdSe33_"+str(subtraj)+"_SP_Dyn_"+str(istate)+".png", dpi=300)
    #start_time += subtraj_increment





pool = mp.Pool(nsubtrajs)
pool.map( myfunc, list(range(nsubtrajs)) )
pool.close()
pool.join()

