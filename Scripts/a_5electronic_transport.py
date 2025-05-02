import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import time
import pickle
from scipy.interpolate import interpn
from scipy.optimize import curve_fit
from numpy import linalg as LA
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})


def get_trjs(file_name):
    file1 = open(file_name, 'r')
    trjs = 0
    time = []
    line = file1.readline()
    while line != '':
        if 'ITEM: TIMESTEP' in line:
            line = file1.readline()
            time.append(int(line))
            trjs += 1
        line = file1.readline()
    times = np.array(time)
    return(trjs,times)

def get_ovito(file_name,trjs):
    file1 = open(file_name, 'r')
    for i in range(0,3):
        file1.readline()
    n = int(file1.readline())    
    file1.readline()
    box = np.array([0,0,0.0])
    for i in range(3):
        line = file1.readline()
        line = line.split()
        box[i] = float(line[1])-float(line[0])
    file1.seek(0)
    ovito = np.zeros((trjs,n,12), dtype = "float")
    for i in range(trjs):
        for j in range (0, 9):
            file1.readline()
        for j in range (0,n):
            line =  file1.readline()
            line = line.split()
            for k in range(0,12): 
                ovito[i,j,k]=float(line[k])
    file1.close()
    return(ovito,box)
def get_nbb(ovito):
    nreal = ovito.shape[0]
    nbb = 0
    for j in range(nreal):
        if ovito[j,1] == 1: #if atom is an ellipsoid
            nbb += 1
    return(nbb)

def get_ellipsoids(ovito):
    nreal = ovito.shape[0]
    nbb = 0
    for j in range(nreal):
        if ovito[j,1] == 1: #if atom is an ellipsoid
            nbb += 1
    ellipsoids = np.zeros((nbb,12),dtype=float)
    k = 0
    for j in range(nreal):
        if ovito[j,1] == 1: #if atom is an ellipsoid
            ellipsoids[k,:] = ovito[j,:12]
            k += 1
    ellipsoids = np.delete(ellipsoids,(0,1,5,6,7,8,9,10,11),axis = 1)
    return(ellipsoids)

def get_center(eigvec,rBB,box):
    n_max = eigvec.argmax()
    rBB_temp = rBB-rBB[n_max] #shift center of box
    rBB_temp = rBB_temp-box*np.around((rBB_temp)/box) #periodic boundaries
    center = eigvec[:,np.newaxis] * rBB_temp  #'mass' of each particle
    center = np.sum(center,axis=0) #center of mass
    center += rBB[n_max] #in original box coordinates
    return(center)

def get_Pis(k_mat, dopants, threshold, rounds,output,eigvals,f_mag,r_mat,PisFile,f):
    n = np.shape(k_mat)[0]
    if PisFile == '':
        pop0 = get_pop0(k_mat,dopants,eigvals)
    else:
        pop0 = np.zeros(n,dtype=float)
        fileo = open(PisFile,'rb')
        pop_old = pickle.load(fileo)
        fileo.close()
        n_old = np.shape(pop_old)[0]
        pop0[-n_old:] = pop_old
    pop1 = np.array(pop0) #new guess
    norm = np.sum(k_mat,axis = 1)
    denom = np.zeros((n,n),dtype = float)
    for i in range(n):
        denom[i] = (k_mat[i,:]-k_mat[:,i])        
    d_pop = 1 #max change in population
    mob0 = 0
    roundj = 0 #current round
    while d_pop > threshold and roundj < rounds:
        roundj += 1
        for i in range(n): #from Bishop 2001
            pop1[i] = np.sum(pop1*k_mat[:,i])/norm[i]/(1-(np.sum(denom[i]*pop1)/norm[i]))
        pop1 = pop1/np.sum(pop1)*dopants #enforce a constant number of dopants
        if roundj%output == 0:
            change = np.abs(pop1-pop0)
            d_pop = np.amax(change)
            pop0 = np.array(pop1)
            print('Iter {}, delta = {:.2f}'.format(roundj,np.log10(d_pop)))
            mobility = get_mobility(pop0, k_mat, f_mag, r_mat,f)
            d_mob = mobility-mob0
            mob0 = mobility
    return(pop1,mobility,d_pop/output,d_mob)

def get_pop0(k_mat,dopants,eigvals):
    kbT = 0.02585 #eV - 300K
    exp = np.exp(-1*-1*eigvals/kbT) #-1 to account for holes
    z = np.sum(exp)
    pop0 = exp/z*dopants
    return(pop0)

def get_mobility(Pis, k_mat, f_mag, r_mat,f):
    matrix = k_mat*Pis[:,np.newaxis]*(1-Pis[np.newaxis,:])
    velocity = np.sum(matrix[:,:,np.newaxis]*r_mat,axis = (0,1))
    mobility = velocity[f]/f_mag
    return(mobility)
