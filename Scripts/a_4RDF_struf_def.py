import numpy as np
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
import math
import pickle
from scipy import integrate

def get_ovito(file_name,trjs):
    file1 = open(file_name, 'r')
    for i in range(0,3):
        file1.readline()
    n = int(file1.readline())
    file1.readline()
    box = np.array([0.0,0.0,0.0])
    for i in range(3):
        line = file1.readline()
        line = line.split()
        box[i] = float(line[1])-float(line[0])
    for i in range(n+3):
        line = file1.readline()
    steps = int(line)
    file1.seek(0)
    ovito = np.zeros((trjs,n,4), dtype = "float")
    for i in range(trjs):
        for j in range (0, 9):
            file1.readline()
        for j in range (0,n):
            line =  file1.readline()
            line = line.split()
            for k in range(0,4): 
                ovito[i,j,k]=float(line[k+1])
    file1.close()
    return(ovito,box,n,steps)

def RDF_atoms(ovito, atom1, atom2):
    atom1del = []
    atom2del = []
    for i in range (ovito.shape[1]):
        if ovito[0,i,0] != atom1:
            atom1del.append(i)
        if ovito[0,i,0] != atom2:
            atom2del.append(i)
    atom1s = np.delete(ovito, atom1del, 1)
    atom1s = np.delete(atom1s, 0, 2) #removes atom id
    atom2s = np.delete(ovito, atom2del, 1)
    atom2s = np.delete(atom2s, 0, 2) #removes atom id
    return(atom1s, atom2s)

def pairs(atom1s, atom2s, d_r, nbins,box): #creates hist from r (minimal image)
    trjs = atom1s.shape[0]
    n = atom1s.shape[1]
    rij_mag = np.zeros((trjs,n,n),dtype = float)
    for i in range(trjs): #Breaking this up speeds up the code substantially.
        atom1s[i,:,:] = atom1s[i,:,:] / box
        atom2s[i,:,:] = atom2s[i,:,:] / box
        rij = atom1s[i,:,np.newaxis,:]-atom2s[i,np.newaxis,:,:] #all distance vectors
        rij = rij - np.rint(rij) #minimal image  ##This step is slow
        rij = rij*box
        rij_mag[i] = np.sqrt(np.sum(rij**2,axis = -1)) #sum xyz  ##This step is slow
        if i%5 == 0:
            print('Trj {} of {} done'.format(i+1,trjs))
    hist = np.zeros((trjs,nbins),dtype = float)
    for i in range(trjs):
        hist[i] , edges = np.histogram(rij_mag[i],bins=nbins,range=(0.0,(d_r*nbins)))
    hist[:,0] = 0 #avoids spike at 0 if atom1 = atom2
    return(hist,edges)

def norm_hist(all_hist, steps, box, d_r, nk, n1, n2): #normalizes the histogram
    hist = np.sum(all_hist, axis = 0)
    hist = hist/(n1*steps)/(4/3*np.pi)/(n2/(box[0]*box[1]*box[2]))
    for i in range (0,nk):
        hist[i] = hist[i]/(((i+1)*d_r)**3-(i*d_r)**3)
    return(hist)

def plot_rdf(hist,atom1,atom2,xmax,start,stop,d_r): #plots the pair dist. function
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"
    fig, ax = plt.subplots(figsize=(4,4))
    #alphas = np.flip(np.arange(0.2,1.0,0.1))
    colors = np.array([(0.1,0.9,0.1)],dtype = float)
    ax.plot(hist[0], hist[1], label='{}:{}'.format(atom1,atom2), color=colors[0])#, alpha = alphas[i]) 
    ax.set_ylabel('RDF', fontsize='x-large')
    ax.set_xlabel('Distance (sigma)', fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.axhline(y=1, color="gray",linewidth=0.1)
    ax.legend(title = "Atoms", ncol=2)
    #ax.text(12.2,47, '$\infty$',zorder=10,fontsize=16)
    ax.set(xlim=(0,30), ylim=(0,10))
    fig.tight_layout()
    plt.savefig('Ave_RDF{}-{}_dr{}.png'.format(start,stop,d_r), dpi=200)
    return


def save_data(pickle_name,all_hist):
    rdffile = open(pickle_name,'wb')
    pickle.dump(all_hist,rdffile)
    rdffile.close()
    return

def load_data(pickle_name):
    rdffile = open(pickle_name,'rb')
    all_hist = pickle.load(rdffile)
    rdffile.close()
    return(all_hist)

def rdf2struf(rdf,box,n1,kmax):
    density = n1/(box[0]*box[1]*box[2])
    boxmin = np.amin(box)/(2*math.pi)/2
    integral = np.zeros((2,int(kmax*boxmin)),dtype = float)
    for j in range(len(integral[0])):
        integral[0,j] = (j+1)/boxmin
    for j in range(len(integral[0])): #k values
        k = integral[0,j]
        fr = rdf[0]**2*(rdf[1]-1)*np.sin(k*rdf[0])/(k*rdf[0])
        integral[1,j] = integrate.simpson(fr,rdf[0])
    integral[1] = 1 + 4*math.pi*density*integral[1]
    return(integral)

def plot_struf(struf,start,stop):
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(struf[0],struf[1])
    ax.set_ylabel('Structure Factor', fontsize='x-large')
    ax.set_xlabel('q ($\sigma^{-1}$)', fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.axhline(y=0, color="gray",linewidth=0.1)
    ax.axhline(y=1, color="gray",linewidth=0.1)
    #ax.legend(title = "Atoms", ncol=2)
    #ax.text(1,1,'kmin={}\nkwidth={}'.format(kmin,kwidth),fontsize=8)
    #ax.set(xlim=(0,8), ylim=(0,2))
    fig.tight_layout()
    plt.savefig('Ave_Struf{}-{}.png'.format(start,stop),dpi = 200)
    return()

def plot_lnstruf(struf,start,stop):
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"
    fig, ax = plt.subplots(figsize=(3,3))
    ax.loglog(struf[0],struf[1], marker = '2', color='slateblue')
    #ax.set_yscale("log", base=2)
    ax.set_ylabel('Structure Factor', fontsize='x-large')
    ax.set_xlabel('q ($\sigma^{-1}$)', fontsize='x-large')
    ax.tick_params(labelsize='large')
    #ax.axhline(y=0, color="gray",linewidth=0.1)
    ax.axhline(y=1, color="gray",linewidth=0.1)
    #ax.legend(title = "Atoms", ncol=2)
    #ax.text(1,1,'kmin={}\nkwidth={}'.format(kmin,kwidth),fontsize=8)
    #ax.set(xlim=(0,8), ylim=(0,2))
    fig.tight_layout()
    plt.savefig('Ave_loglogStruf{}-{}New.png'.format(start,stop),dpi = 200)
    return()

