import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import sys
import os
import time
from a_4ovito import *
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

#Changing inputs
files = [i for i in os.listdir(".")]
if 'ovito.trj' not in files:
    print('No ovito.trj')
    exit()
filename = 'ovito.trj'
chainlength = 64 #monomers
atoms_pm = 5 #atoms per monomer in ovito

trjs,steps = get_trjs(filename)
if 'savefile_run.out' not in files:
    time_step = 0.002 #tau
    print('No savefile_run.out, assuming a {} tau timestep'.format(time_step))
else:
    time_step = get_timestep('savefile_run.out')
    print('Using a {} tau timestep based on savefile_run.out'.format(time_step))
times = steps*time_step #converts to units of tau
trj_int = times[1] - times[0]

#Get ovito.trj data
ovito_raw = get_ovito(filename, trjs)
ellipsoids_raw = get_ellipsoids(ovito_raw)
anions_raw = get_anions(ovito_raw)
ovito, chains = break_chains(ovito_raw, chainlength, atoms_pm)
ellipsoids, chains = break_chains(ellipsoids_raw, chainlength, 1) #one ellipsoid per chain

"""pi_orient, bond_orient, dihedrals = calc_dihedrals(ellipsoids)
#orient = [trj, chain, n (or n-1), 3(xyz)], dihedrals = [trj, chains, n-1]
#pi and bond are normalized
#dihedrals are from [-pi,pi]

dihedrals_raw = chains_to_raw(dihedrals)
pi_orient_raw = chains_to_raw(pi_orient)
bond_orient_raw = chains_to_raw(bond_orient)
#orient_raw = [trj, n_total, 3], dihedrals_raw = [trj, n_total]
"""
#############################################################
# Output Options: ovito (_raw), ellipsoids (_raw), times    #
# ovito & ellipsoids matrix is: [trj, chain, atom, 12]      #
# anions_raw [trj, atom, 12]                                #
# dihedrals (_raw), pi_orient (_raw), bond_orient (_raw)    #
# [trj, chain, n n-1, 1 xyz] or _raw [trj, chain*n, 1 xyz]  #
#############################################################

starts = [0]
stops = [51]
Temps = [1]
out = []

for i in range (len(starts)):
    start = starts[i] #1st trj to include in analysis
    stop = stops[i] #last trj to include in analysis
    temp = Temps[i]
    out.append('ForTRJs {}-{}\n'.format(start,stop))
    """
    plot_dihedrals(dihedrals_raw[start:stop,:],start,stop,temp)

    time_corl_max = 10 #end of the time_corl graph in units of trj_int
    dihed_time_corl = calc_dihedral_time_corl(dihedrals_raw[start:stop,:],time_corl_max,trj_int)
    plot_dihedral_time_corl(dihed_time_corl,start,stop) #dihedxy_ti.dihedxy_ti+j vs. time interval(j)

    bond_corl_max, corl_fit_i,corl_fit_f = 60, 1, 10 #max bond spacing analyzed, Max N-2 
    Lp,B_max = plot_bondcorrelation(ellipsoids[start:stop],start,stop,bond_corl_max,corl_fit_i,corl_fit_f) #bondi.bondi+1 vs. bonds
    out.append('L_p(bonds) = {:.4f}\nBonds_Analyzed = {:.0f}'.format(Lp,B_max))
    """
    ### Axis not edited to work ###
    #axis_corl_max, axis_op_max, d_r = 60, 10, 1 #max bond spacing analyzed, Max N-3
    #plot_axis_corl(ellipsoids,start,stop,axis_corl_max) #axis.axis vs. bonds
    #plot_axis_op(ellipsoids,start,stop,axis_op_max, d_r) #axis.axis vs. distance
    
    #Nematic OPs takes 1 second per trj
    """nematic_op_max, d_r = 10, 1 #max distance in sigma analyzed, Max Box/2 (15); histogram spacing
    plot_pinematic_op(pi_orient_raw[start:stop],ellipsoids_raw[start:stop],start,stop,nematic_op_max, d_r) #pi.pi vs. distance
    plot_bondnematic_op(bond_orient_raw[start:stop],ellipsoids_raw[start:stop],start,stop,nematic_op_max, d_r) #b.b vs. distance"""
    """out.append(nemOPs(pi_orient_raw[start:stop],bond_orient_raw[start:stop],start,stop))
    
    #Radius of gyration time decorelation
    Rgs, Rgs_ave, Rgs_std = get_Rgs(ellipsoids[start:stop]) #Rgs = [trj,chain,(rg_total, aspher, acyl, aniso)]
    time_corl_max, corl_fit_i, corl_fit_f = 50, 2, 15 #end of the time_corl graph in units of trj_int
    rg_time_corl = calc_rg_time_corls(Rgs[:],time_corl_max,trj_int) #out = [time_int,chain,(time,rg,aspher,acyl,aniso)]
    Rg_Tcor = plot_rg_time_corl(rg_time_corl,start,stop,corl_fit_i,corl_fit_f)
    out.append('Rg_ave(sig) = {:.2f}\nRg_std(sig) = {:.2f}\nRg_Tcor(tau) = {:.0f}\nAs_ave(sig2) = {:.2f}\nAs_std(sig2) = {:.2f}'.format(Rgs_ave[0],Rgs_std[0],Rg_Tcor,Rgs_ave[1],Rgs_std[1]))  
    
    time_corl_max, corl_fit_i, corl_fit_f = 50, 2, 15 #end of the time_corl graph in units of trj_int
    ends_length_time_corl = get_ends_length_time_corl(ellipsoids[start:stop],time_corl_max,trj_int) #timecorl = [time, (time, chains' e.e)]
    E2El_Tcor = plot_E2El_time_corl(ends_length_time_corl,start,stop,corl_fit_i,corl_fit_f)
    time_corl_max, corl_fit_i, corl_fit_f = 50, 2, 40 #end of the time_corl graph in units of trj_int
    ends_vector_time_corl = get_ends_vector_time_corl(ellipsoids[start:stop],time_corl_max,trj_int) #timecorl = [time, (time, chains' e.e)]
    E2Ev_Tcor = plot_E2Ev_time_corl(ends_vector_time_corl,start,stop,corl_fit_i,corl_fit_f)
    out.append('E2El_Tcor(tau) = {:.0f}\nE2Ev_Tcor(tau) = {:.0f}\n'.format(E2El_Tcor,E2Ev_Tcor))
    """
    #Polymer Diffusion#
    MSD_max, D_fit_i, D_fit_f = 20, 5, 20 #end of the time_corl graph in units of trj_int
    MSD = calc_MeSqDisp(ellipsoids[start:stop], MSD_max, trj_int) #MSD = [time,chains,(time,r2,x2,y2,z2)]
    Diff, Diff0 = plot_MeSqDisp(MSD,start,stop,D_fit_i,D_fit_f)
    out.append('PolyDiff(sig2/ktau) = {:.4f}\n'.format(Diff*1000))
    out.append('PolyDiff0(sig2/ktau) = {:.4f}\n'.format(Diff0*1000))

    #Anion Diffusion#
    MSD_max, D_fit_i, D_fit_f = 20, 2, 10 #end of the time_corl graph in units of trj_int
    MSD = calc_Anion_MeSqDisp(anions_raw[start:stop], MSD_max, trj_int) #MSD = [time,chains,(time,r2,x2,y2,z2)]
    Diff,Diff0 = plot_Anion_MeSqDisp(MSD,start,stop,D_fit_i,D_fit_f)
    out.append('AnionDiff(sig2/ktau) = {:.2f}\n'.format(Diff*1000))
    out.append('AnionDiff0(sig2/ktau) = {:.2f}\n'.format(Diff0*1000))

writeout(out, 'zout_ovito_diffusionOnly.out')

