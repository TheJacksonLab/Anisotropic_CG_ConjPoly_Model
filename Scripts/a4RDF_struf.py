import numpy as np
import matplotlib.pyplot as plt
import math
import os
from a_4RDF_struf_def import *

filename = 'ovito0-51.trj' #input("Enter file name: ")
atom1 = 1
atom2 = 1
trj_initial = [0] #int(input("What is the first time step? (>=0): "))
trj_final = [51] #int(input("What is the last time step? : "))
trjs = trj_final[0]-trj_initial[0]

ovito, box, n, steps = get_ovito(filename, trj_final[0]) #ovito = atom type, x, y, z (sorted by atom id)
timestep = 0.002
print('Assuming a {} tau timestep'.format(timestep),flush=True)

atom1s, atom2s = RDF_atoms(ovito[trj_initial[0]:trj_final[0]], atom1, atom2) #two arrays of atom1s & 2s (
n1 = atom1s.shape[1]
n2 = atom2s.shape[1]

##### Calculate rij Histograms #####
d_r = 0.01 # the width of each bin in sigma
length = np.amin(box)/2
nbins = int(np.floor(length/d_r))

pickle_name = 'all_hist{}-{}_dr{}.pkl'.format(trj_initial[0],trj_final[0],d_r)
all_hist, edges = pairs(atom1s,atom2s, d_r, nbins, box) #this is the slow step in the code
save_data(pickle_name,all_hist)
all_hist = load_data(pickle_name)

## Time Average RDF ##
ave_rdf = np.zeros((2,nbins),dtype = float)
ave_rdf[0] = np.arange(d_r/2,d_r*nbins+d_r/2,d_r)
ave_rdf[1] = norm_hist(all_hist, trjs, box, d_r, nbins, n1, n2)
plot_rdf(ave_rdf,atom1,atom2,d_r*nbins,trj_initial[0],trj_final[0],d_r)

## Time Average Structure Factor ##
kmax = 2*math.pi*1.5 #This keeps the x-axis max = 10 
struf = rdf2struf(ave_rdf,box,n1,kmax)
plot_struf(struf,trj_initial[0],trj_final[0])
plot_lnstruf(struf,trj_initial[0],trj_final[0])
