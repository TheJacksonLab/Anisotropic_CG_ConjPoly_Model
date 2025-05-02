import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import math
import sys
import time
from a_5electronic_transport import *
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

tmax = 50
ovitofile = 'ovito_100_150.trj'
#trjs, steps  = get_trjs(ovitofile)
ovito_raw, box = get_ovito(ovitofile, tmax)
n = get_nbb(ovito_raw[0])
rBB_all = np.zeros((tmax,n,3),dtype = float)
for t in range(tmax):
    rBB_all[t] = get_ellipsoids(ovito_raw[t]) #[nbb,3(x,y,z)]

ham_all = np.load("H100-150.npy")
eigvals_all = np.load("eigenvalues100-150.npy")
eigvecs_all = np.load("eigenvectors100-150.npy")

print("Pickles Loaded")
if n != ham_all.shape[1]:
    print("Number of beads in simulation != Dimension of hamiltonian!")
    exit

f_mag = 1e-5  #electric field strength eV/sigma
f_vecs = f_mag*np.array([[1,0,0],[0,1,0],[0,0,1]])

mobility = np.zeros((tmax,np.shape(f_vecs)[0]), dtype = float)
d_pop = np.zeros((tmax,np.shape(f_vecs)[0]), dtype = float)
d_mob = np.zeros((tmax,np.shape(f_vecs)[0]), dtype = float)

for t in range(tmax):
    print('Trj = {}'.format(t), flush=True)
    ham = ham_all[t]
    eigvals = eigvals_all[t]
    eigvecs = eigvecs_all[t]
    rBB = rBB_all[t]
    
    ## Calculate kij Matrix ##
    lam_coup = 0.45 #eV
    g_coup = 0.005 #unitless
    kbT = 0.02585 #eV, 300K
    hbar = 6.58212e-16 #eV*s
    ham *= (np.ones((n,n),dtype = int) - np.identity(n)) #set diagonal = 0
    s_mat = ((eigvecs.T)**2)@ham**2@(eigvecs**2) #eV**2
    s_mat *= (np.ones((n,n),dtype = int) - np.identity(n)) #set diagonal = 0 CHECKED
    
    a_temp = np.array([np.sum(eigvecs**4,axis=0)]) #IPR per state CHECKED
    lam_mat = lam_coup*(a_temp.T + a_temp)
    
    #print(np.amax((s_mat)*g_coup**2/lam_mat**2)) #Confirm Troisi 2014 Eq 22 is True
    #print(np.average((s_mat)*g_coup**2/lam_mat**2))
    
    e_temp = np.array([eigvals])
    e_mat = -1*(e_temp - e_temp.T) #CHECKED
    # '-1*' b/c I'm lookin at a hole in basis of electron HOMOs

    cen_vec = np.zeros((n,3),dtype = float)
    for i in range(n):
        cen_vec[i] = get_center(eigvecs[:,i]**2,rBB,box) #centers of each MO CHECKED
    r_mat = cen_vec[np.newaxis,:,:]-cen_vec[:,np.newaxis,:]
    r_mat -= box*np.around(r_mat/box) #r_mat(n[i],n[j],3(xyz)) CHECKED
    for f in range(np.shape(f_vecs)[0]):
        f_mat = np.sum(r_mat*f_vecs[f],axis=2) #CHECKED on homo-1 to homo-2
        
        #This is the hopping rate matrix between energy levels, just a lot of elementwise operations
        #on the matrices previously set up. k_mat = [from state, to state]
        k_mat = g_coup**2*s_mat/hbar*np.sqrt(np.pi/kbT/lam_mat)*np.exp(-((lam_mat+e_mat-f_mat)**2)/(4*kbT*lam_mat))
        #holes flow up field
        #print(k_mat[-3:,-3:])
        
        rounds = 200000 #max rounds
        dopants = 10
        threshold_true = 1e-8 #max allowable change per step
        output = 1000
        threshold = threshold_true*output #max allowable change per output
        ns = 1000
        PisFile = '' #pkl file of another run if continuing equilibration
        #starttime = time.time()
        Pis,mobility[t,f],d_pop[t,f],d_mob[t,f] = get_Pis(k_mat[-ns:,-ns:], dopants, threshold, rounds, output,eigvals[-ns:],f_mag, r_mat[-ns:,-ns:],PisFile,f)
        #print(time.time()-starttime)
        dim = ['x','y','z']
        fileo = open('Pis_t{}{}.pkl'.format(t,dim[f]),'wb')
        pickle.dump(Pis, fileo)
        fileo.close()
        print(mobility[t,f],flush=True)

f1 = open('zout_mobility_F{:.2f}_D{}.csv'.format(np.log10(f_mag),dopants),'w')
f1.write('Mob_ave,Mob_sig,mob_x,mob_y,mob_z,sig_x,sig_y,sig_z\n')
f1.write('{},'.format(np.average(mobility)))
f1.write('{},'.format(np.std(mobility)))
for f in range(3):
    f1.write('{},'.format(np.average(mobility[:,f])))
for f in range(2):
    f1.write('{},'.format(np.std(mobility[:,f])))
f1.write('{}\n\n'.format(np.std(mobility[:,f+1])))
f1.close()

f1 = open('zout_fullmobility_F{:.2f}_D{}.csv'.format(np.log10(f_mag),dopants),'w')
f1.write('TRJ,Mob_ave,Mobx,Moby,Mobz,d_mobx,d_moby,d_mobz,log(d_popx),log(d_popy),log(d_popz)\n')
for t in range(np.shape(mobility)[0]):
    f1.write('{},'.format(t))
    f1.write('{},'.format(np.average(mobility[t])))
    for f in range(np.shape(mobility)[1]):
        f1.write('{},'.format(mobility[t,f]))
    for f in range(np.shape(mobility)[1]):
        f1.write('{},'.format(d_mob[t,f]/mobility[t,f]))
    for f in range(np.shape(mobility)[1]-1):
        f1.write('{},'.format(np.log10(d_pop[t,f])))
    f1.write('{}'.format(np.log10(d_pop[t,f+1])))
    f1.write('\n')
f1.close()
