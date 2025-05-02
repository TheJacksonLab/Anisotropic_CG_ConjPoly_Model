import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import math
import sys
import time
from a_4bundles_def import *
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

#Changing inputs
ovitofile = 'ovito.trj'
trjs, steps  = get_trjs(ovitofile) #trjs is an int
time_step = 0.002
times = steps*time_step #converts to units of tau, times is a vector
trj_int = times[1] - times[0]
chainlength = 64 #monomers per chain
atoms_pm = 5 #atoms per monomer in ovito

#Get ovito.trj data
ovito_raw, box = get_ovito(ovitofile, trjs)
ellipsoids_raw = get_ellipsoids(ovito_raw)
ellipsoids, chains = break_chains(ellipsoids_raw, chainlength, 1) #one ellipsoid per monomer

Aij = np.zeros((trjs,chainlength*chains,chainlength*chains), dtype=bool)
A_cut = 1.7
exclude = 5
maxblock = 3
out = ['For A_cut = {}, exclude = {}'.format(A_cut,exclude)]

start, stop = 10,51
trjs = stop-start
bundlename = 'Bundles{}-{}_{}cut_ex{}.pkl'.format(start,stop,A_cut,exclude)
d_r, nbins = 1, 30
RDFname = 'RDFs{}-{}.pkl'.format(start,stop)

### Calculate Bundles ### Comment out for reruns ###
for i in range(trjs):
    if i%10 == 0:
        print('Aij at trj {} of {}'.format(i,trjs))
    Aij[i,:,:] = get_Aij(chains,chainlength,ellipsoids_raw[start+i,:,2:5],box,A_cut)

bundles= np.zeros((trjs,Aij.shape[1]),dtype = float)
for i in range(trjs):
    bundles[i] = get_bundles(Aij[i],exclude,maxblock,chainlength)
    if i%10 == 0:
        print('Bundles at trj {} of {}'.format(i,trjs))

efile = open(bundlename, 'wb')
pickle.dump(bundles,efile)
efile.close()

## Correlation Length via g_inter vs. g_intra ##
RDFs = np.zeros((trjs,4,nbins),dtype = float)
RDFs[:,0,:] = np.arange(0,nbins,d_r)+0.5
for i in range(trjs):
    RDFs[i,1:3] = get_intrainterRDFs(ellipsoids_raw[start+i,:,2:5],Aij[i],exclude,maxblock,chainlength,box,d_r,nbins)
    if i%10==0: print('RDF {}/{} complete'.format(i+1,trjs))
RDFs[:,3,:] = RDFs[:,1,:]+RDFs[:,2,:]

efile = open(RDFname, 'wb')
pickle.dump(RDFs,efile)
efile.close()
#### End comment out #### End comment out ####

efile = open(bundlename, 'rb')
bundles = pickle.load(efile)
efile.close()

efile = open(RDFname, 'rb')
RDFs = pickle.load(efile)
efile.close()

## Write Bundled TRJ ##
fileout = 'ovitobundles{}-{}_{}cut_exclude{}.trj'.format(start,stop,A_cut,exclude)
write_trjwithbundles(ovitofile,bundles,fileout,start)

## RDF and Corl Length ##
CorlLeng = np.zeros(trjs,dtype = float)
for i in range(trjs):
    CorlLeng[i] = get_CorlLeng(RDFs[i]) #linear interpolation between points

Ave_RDF = np.average(RDFs,axis = 0)
plot_RDFs(Ave_RDF,d_r*nbins,start,stop,np.average(CorlLeng))
out.append('Corl_Length(sigma) = {:.0f}'.format(np.average(CorlLeng)))

## Histograms ##
bundle_ave, bundle_std = bundle_histogram(bundles[0:trjs],trj_int,start,stop)
out.append('Bundle_Ave(monomer_weighted) = {:.2f}\nBundle_Std(monomer_weighted) = {:.2f}'.format(bundle_ave,bundle_std))

bundle_histogram_time(bundles[0:trjs],trj_int,start,stop)

## Time Correlation ##
corl_max, corl_fit_i, corl_fit_f = 40, 1, 30
ave_time_corl = bundle_time_corl(bundles[0:trjs],trj_int,corl_max)
bundle_Tcor = plot_bundle_time_corl(ave_time_corl,start,stop,corl_fit_i,corl_fit_f) 
out.append('Bundle_Tcor(tau) = {:.0f}'.format(bundle_Tcor))

writeout(out,'{}cut_exclude{}.out'.format(A_cut,exclude))
