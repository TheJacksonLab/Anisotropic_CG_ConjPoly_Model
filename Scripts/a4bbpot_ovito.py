import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import math
import sys
import time
from a_4bbpot_ovito import *
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

#Changing inputs
start, stop = 0, 51 
ovitofile = 'ovito{}-{}.trj'.format(start,stop)
nbbpotsfile = 'nbbpots{}-{}.out'.format(start,stop)
trjs, steps  = get_trjs(ovitofile) #trjs is an int
time_step = 0.002
times = steps*time_step #converts to units of tau, times is a vector
trj_int = times[1] - times[0]
chainlength = 64 #monomers per chain
atoms_pm = 5 #atoms per monomer in ovito

#constants
rij_decay = 1 #/sigma - from Savoie 2021
rij_0 = 0.75 #sigma - based on minimum rij
J_inter = 0.1 #eV -guessed from Troisi 2007 and Venkataraan 2012, 5x higher than Savoie 2021
e_diag = 0  #eV to kcal/mol, 2eV in 2016 paper "Tight Binding... polythiophenes" 
e_off = -1  #eV, from 2016 paper "Tight Binding... polythiophenes"
#note: H has -e_off on the off diagonal elements)

#Get ovito.trj data
ovito_raw, box = get_ovito(ovitofile, trjs)
ellipsoids_raw = get_ellipsoids(ovito_raw)
ovito, chains = break_chains(ovito_raw, chainlength, atoms_pm)
ellipsoids, chains = break_chains(ellipsoids_raw, chainlength, 1) #one ellipsoid per chain

pi_orient, bond_orient, dihedrals = calc_dihedrals(ellipsoids)
#orient = [trj, chain, n (or n-1), 3(xyz)], dihedrals = [trj, chains, n-1]
#pi and bond are normalized
#dihedrals are from [-pi,pi]

dihedrals_raw = chains_to_raw(dihedrals)
pi_orient_raw = chains_to_raw(pi_orient)
bond_orient_raw = chains_to_raw(bond_orient)
#orient_raw = [trj, n_total, 3], dihedrals_raw = [trj, n_total]

#Get bbpots
bbpots = read_bbpots(nbbpotsfile, trjs)

Cz = 1 #int(sys.argv[1]) #couplings
Dz = 1 #int(sys.argv[2]) #dihedrals
Ez = 1 #int(sys.argv[3]) #electrostatics
monos = chainlength*chains
eigvecs = np.zeros((trjs,monos,monos), dtype=float)
eigvals = np.zeros((trjs,monos), dtype = float)
h_off = np.zeros((trjs,monos-1), dtype =float)
wij = np.zeros((trjs,monos,monos), dtype = float)
H = np.zeros((trjs,monos,monos), dtype=float)

for i in range(trjs):
    if i%1 == 0:
        print('At trj {} of {}'.format(i+1,trjs))
    eigvecs[i,:,:], eigvals[i,:], h_off[i,:], wij[i,:], H[i,:,:] = get_eigs(i,chains,chainlength, e_diag, e_off,bbpots[i], dihedrals_raw[i,:],ellipsoids_raw[i,:,2:5],pi_orient_raw[i,:,:],J_inter,rij_decay,rij_0, Cz,Dz,Ez,box)

efile = open('eigvals{}-{}.pkl'.format(start,stop), 'wb')
pickle.dump(eigvals,efile)
efile.close()
efile = open('eigvecs{}-{}.pkl'.format(start,stop), 'wb')
pickle.dump(eigvecs,efile)
efile.close()
efile = open('ham{}-{}.pkl'.format(start,stop), 'wb')
pickle.dump(H,efile)
efile.close()

#############################################################
# Output Options: eigvals, eigvecs, bbpot, ellipsoids, dihedrals, h_off, wij(coupling)#
#############################################################
"""
fileo = open('ham{}-{}.pkl'.format(start,stop), 'rb')
ham = pickle.load(fileo)
fileo.close()
fileo = open('eigvals{}-{}.pkl'.format(start,stop), 'rb')
eigvals = pickle.load(fileo)
fileo.close()
fileo = open('eigvecs{}-{}.pkl'.format(start,stop), 'rb')
eigvecs = pickle.load(fileo)
fileo.close()
print('Pickle Files Loaded')
"""
regions = 1 #analysis regions
starts = [0]
stops = [51]
Temps = [1]
out = []

for j in range (0,regions):
    trjs = stops[j]-starts[j]
    start = starts[j]
    stop = stops[j]
    out.append('ForTRJs {}-{}'.format(start,stop))
    res = 30 #histogram resolution

    #Properties of the system before the Hamiltonian
    out.append(plot_deltabbpot(bbpots[start:stop],res,start,stop))
    out.append(bbpot_histogram(bbpots[start:stop],res,start,stop))

    ###bbpot_histogram2(bbpots[start:stop],res,start,stop)
    out.append(print_h_off_stats(h_off[start:stop,:],start,stop))
    neighbors = 2
    plot_wij(wij,res,start,stop,neighbors)
        
writeout(out, 'zout_bbpot_ovito.out')

#Save a new trj with the bbpot and HOMO
MO0 = chains*chainlength-10 #chains*chainlength-1 = just HOMO
MOf = chains*chainlength #chains*chainlength means HOMO is the last orbital plotted
write_trjwithMOs(ovitofile,eigvecs[:,:,MO0:MOf],MO0,MOf,bbpots)
