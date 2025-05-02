import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import math
import os
import sys
import time
from a_4bbpot_percolate import *
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

#make_ham()
trj_filename = 'ovito0-51.trj'
fileoname = ['ham0-51.pkl'] #i for i in os.listdir(".") if ".pkl" in i and 'ham' in i]
fileo = open(fileoname[0], 'rb')
ham = pickle.load(fileo)
##ham's dimension was decreased
#ham = ham[np.newaxis,:,:]
###############################
fileo.close()
print('Pickled Ham Imported')


#############################################################
# Output Options: eigvals, eigvecs, bbpot, ellipsoids, dihedrals, h_off, wij(coupling)#
#############################################################

regions = 1 #analysis regions
starts = [0]
stops = [51]
out = []

for j in range (0,regions):
    trjs = stops[j]-starts[j]
    start = starts[j]
    stop = stops[j]
    out.append('ForTRJs {}-{}'.format(start,stop))

    threshold = 0.01 #eV
    netIDs = np.zeros((trjs,ham.shape[1]),int)
    networks = np.zeros((trjs),int)
    for i in range(trjs):
        netIDs[i,:], networks[i] = get_netIDs(ham[i],threshold)
    max_nets = np.amax(networks)
    net_sizes = np.zeros((trjs,max_nets),dtype = int) #Monos per network, sorted; index is Network number
    for i in range(trjs):
        netIDs[i,:],net_sizes[i,:(networks[i])] = sort_netIDs(netIDs[i,:],networks[i])

write_trjwithnetworks(trj_filename, netIDs, threshold, np.amax(networks))
#print('{} Networks'.format(networks))
writeout(out, 'zout_percolate.out')
