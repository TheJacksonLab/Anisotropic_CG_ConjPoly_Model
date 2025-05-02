import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import time
from scipy.interpolate import interpn
from scipy.optimize import curve_fit
from numpy import linalg as LA
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

def make_ham():
    ham = np.zeros((6,6),int)
    net1 = [0,5]
    net2 = [3,4]
    net3 = [4,2]
    for i in range(len(net1)):
        for j in range(len(net1)):
            ham[net1[i],net1[j]] = 1
    for i in range(len(net2)):
        for j in range(len(net2)):
            ham[net2[i],net2[j]] = 1
    for i in range(len(net3)):
        for j in range(len(net3)):
            ham[net3[i],net3[j]] = 1
    for i in range(ham.shape[0]):
        ham[i,i] = 0
    print(ham)
    return(ham)

def get_netIDs(ham,threshold):
    nmax = ham.shape[0]
    netIDs = np.zeros(nmax,int)
    network = -1
    unassigned = np.ones(nmax,bool)
    innet = np.zeros(nmax,bool)
    adjmat = (ham > threshold)
    n=0
    while n < nmax:
        if unassigned[n] == True:
            network += 1 #start a new network
            innet[n] = True
        while np.sum(innet) != 0:
            i = n
            while i < nmax:
                if innet[i] == True:
                    netIDs[i] = network
                    unassigned[i] = False 
                    innet = (innet | adjmat[i]) & unassigned
                    #i = n-1 #<n have been assigned, but skips can occur
                i += 1
                
        n += 1
    network += 1
    return(netIDs,network)

def sort_netIDs(oldIDs,networks):
    newIDs = np.zeros_like(oldIDs)
    n = oldIDs.shape[0]
    convert = np.zeros((networks,3),dtype=int) #old ID, temp_monos, new ID
    convert[:,0] = np.arange(networks)
    for i in range(n): #make histogram in 
        net = oldIDs[i]
        convert[net,1] += 1
    for i in range(networks): #generate newIDs
        size_max = int(np.amax(convert[:,1]))
        j = 0
        while convert[j,1] != size_max:
            j += 1
        convert[j,2] = i
        convert[j,1] = 0
    for i in range(n): #assign new IDs
        newIDs[i] = convert[oldIDs[i],2]
    for i in range(n): #make new histogram
        net = newIDs[i]
        convert[net,1] += 1
    return(newIDs,convert[:,1])

def write_trjwithnetworks(trj_filename, networkIDs, threshold, maxnetworks):
    f_in = open(trj_filename, 'r')
    f_out = open('ovitoNetwork_{}meV_max{}nets.trj'.format(threshold,maxnetworks),'w')
    trjs = networkIDs.shape[0]
    n = networkIDs.shape[1]
    atomEP = [-1,-1,-1,-1] #all non-BB beads are in network -1
    for i in range(trjs):
        ## TRJ Header ##
        lines = ''
        for j in range(8):
            line = f_in.readline()
            lines += line
        f_out.write(lines)
        line = f_in.readline()
        splitline = line.split()
        splitline.append('network')
        line = ''
        for k in range(len(splitline)):
            line += '{} '.format(splitline[k])
        line += '\n'
        f_out.write(line)
        
        ## TRJ beads ##
        for j in range(n):  #Note: the 1st atom must be a bb bead
            line = f_in.readline()
            splitline = line.split()
            splitline.append(networkIDs[i,j])
            line = ''
            for k in range(len(splitline)):
                line += '{} '.format(splitline[k])
            line += '\n'
            f_out.write(line)
            for k in range(4): #4 sidechain beads
                line = f_in.readline()
                splitline = line.split()
                splitline.append('{}'.format(atomEP[k]))
                line = ''
                for l in range(len(splitline)):
                    line += '{} '.format(splitline[l])
                line += '\n'
                f_out.write(line)
    f_in.close
    f_out.close
    return
"""
def get_n(file_name):
    f = open(file_name, 'r')
    for i in range(0,3):
        f.readline()
    n = int(f.readline())    
    f.close
    return(n)

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

def read_bbpots(file_name,trjs):
    f = open(file_name, 'r')
    for i in range(0,3):
        f.readline()
    n = int(f.readline())    
    f.seek(0)
    bbpots = np.zeros((trjs,n), dtype = "float")
    for i in range(trjs):
        for j in range(0,9):
            f.readline()
        for j in range(0,n):
            line = f.readline()
            line = line.split()
            bbpots[i,j] = float(line[1])
    f.close
    bbpots = bbpots*1000000000*.02585 #converts kT to eV
    return(bbpots)

def get_ovito(file_name,trjs):
    file1 = open(file_name, 'r')
    for i in range(0,3):
        file1.readline()
    n = int(file1.readline())    
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
    return(ovito)

def get_ellipsoids(ovito):
    nreal = ovito.shape[1]
    trjs = ovito.shape[0]
    nbb = 0
    for j in range(nreal):
        if ovito[1,j,1] == 1: #if atom is an ellipsoid
            nbb += 1
    ellipsoids = np.zeros((trjs,nbb,12),dtype=float)
    for i in range(trjs):
        k = 0
        for j in range(nreal):
            if ovito[i,j,1] == 1: #if atom is an ellipsoid
                ellipsoids[i,k,:] = ovito[i,j,:12]
                k += 1
    return(ellipsoids)

def break_chains(data,chainlength,atoms_pm): #
    trjs = data.shape[0]
    n_total = data.shape[1]
    atoms = chainlength*atoms_pm #atoms per chain
    chains = int(n_total/atoms) #number of chains in the box
    if n_total % (atoms) != 0:
        print('The number of atoms does not divide evenly into separate chains')
        exit()
    #else: print('{} chains found'.format(chains))
    data_chains = np.zeros((trjs, chains, atoms, 12), dtype = "float")
    for i in range(chains):
        data_chains[:,i,:,:] = data[:,i*atoms:(i+1)*atoms,:]
    return(data_chains,chains)

def chains_to_raw(data): #data = [trjs,chains,n,x]
    chains = data.shape[1]
    n = data.shape[2]
    if data.ndim == 4: #vector data per monomer
        raw_data = np.zeros(((data.shape[0]),int(chains*n),data.shape[3]),dtype = "float")
        for i in range(chains):
            raw_data[:,(i*n):((i+1)*n),:] = data[:,i,:,:]
    elif data.ndim == 3: #scalar data per monomer
        raw_data = np.zeros(((data.shape[0]),int(chains*n)),dtype = "float")
        for i in range(chains):
            raw_data[:,(i*n):((i+1)*n)] = data[:,i,:]
    return(raw_data)

def rot_mats(quats):
    b = quats[:,:,:,5]
    c = quats[:,:,:,6]
    d = quats[:,:,:,7]
    a = quats[:,:,:,8]
    res = np.zeros((quats.shape[0],quats.shape[1],quats.shape[2],3,3))
    res[:,:,:,0,0] = a*a + b*b - c*c - d*d
    res[:,:,:,0,1] = 2*b*c - 2*a*d
    res[:,:,:,0,2] = 2*b*d + 2*a*c
    res[:,:,:,1,0] = 2*b*c + 2*a*d
    res[:,:,:,1,1] = a*a - b*b + c*c - d*d
    res[:,:,:,1,2] = 2*c*d - 2*a*b
    res[:,:,:,2,0] = 2*b*d - 2*a*c
    res[:,:,:,2,1] = 2*c*d + 2*a*b
    res[:,:,:,2,2] = a*a - b*b - c*c + d*d
    return res

def calc_dihedrals(quats):
    n = quats.shape[2]
    chains = quats.shape[1]
    trjs = quats.shape[0]
    pi_orient = np.zeros((trjs,chains,n,3), dtype = float)
    bond_orient = np.zeros((trjs,chains,n-1,3), dtype = float)
    dihedrals = np.zeros((trjs,chains,n-1),dtype = float)
    rot_maxs = rot_mats(quats[:,:,:,:])
    pi_orient[:,:,:,:] = np.matmul(rot_maxs,[0.0, 1.0, 0.0])
    #pi_orient is the normalized vector of the pi_system
    bond_orient[:,:,:,:] = quats[:,:,1:,2:5] - quats[:,:,:(n-1),2:5]
    bond_norm = LA.norm(bond_orient,axis = -1)
    bond_orient = bond_orient/np.repeat(bond_norm[:,:,:,np.newaxis],3, axis = -1)
    #bond_orient is the normalized vector from a BB atom to the next BB atom
    #https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python  for an explanation
    b0dotb1 = np.sum(pi_orient[:,:,:-1,:]*bond_orient,axis = -1)
    b2dotb1 = np.sum(pi_orient[:,:,1:,:]*bond_orient,axis = -1)
    v = pi_orient[:,:,:-1,:] - np.repeat(b0dotb1[:,:,:,np.newaxis],3,axis = -1)*bond_orient
    w = pi_orient[:,:, 1:,:] - np.repeat(b0dotb1[:,:,:,np.newaxis],3,axis = -1)*bond_orient
    x = np.sum((v*w),axis = -1)
    y1 = np.cross(bond_orient[:,:,:,:],v)
    y = np.sum((y1*w),axis = -1)
    dihedrals = np.arctan2(y,x)
    return(pi_orient, bond_orient, dihedrals)

def get_eigs(trj,chains,chainlength, e_diag, e_off, dihedrals,r,pi,J_inter,rij_decay,rij_0,Cz,Dz,Ez): #Builds Hamiltonian matrix
    n = chains*chainlength
    H = np.zeros((n,n), dtype = float)
    #electronic contribution
    if Ez == 1:
        for i in range(n):
            H[i,i] = e_diag #kT to eV
    h_off = np.zeros(n-1, dtype = float)
    #dihedral contribution (dihedrals only has n-chains elements)
    #cos = math.cos(dihedrals)
    i = 0
    for j in range(chains):
        for k in range(chainlength-1):
            cos = math.cos(dihedrals[(j*(chainlength-1))+k]) #grab the right angle
            if Dz == 1:
                h_off[i] = -e_off*abs(-.01+1.275*cos+.016*cos**2-.87*cos**3-.029*cos**4+.54*cos**5)
            elif  Dz == 0:
                h_off[i] = -e_off
            #coefficients are based on data in Milner 2016 paper
            H[i,i+1] = h_off[i]
            H[i+1,i] = h_off[i]
            i += 1
        i += 1 #skip the coupling between chains
    #coupling contribution
    if Cz == 1:
        rijs = np.zeros((n,n,3), dtype = float)
        rijs = r[np.newaxis,:,:]-r[:,np.newaxis,:]
        rijs_mag = LA.norm(rijs,axis=-1)
        rijs_mag = np.where(rijs_mag == 0, 1, rijs_mag) #put 1s on the diagonal to avoid nans
        rijs = rijs/np.repeat(rijs_mag[:,:,np.newaxis],3,axis=2) #normalized vectors
        wij = np.zeros((n,n), dtype = float) #from Jackson paper
        wij = (np.sum(pi[np.newaxis,:,:]*rijs,axis=2))**2*\
            (np.sum(pi[:,np.newaxis,:]*rijs,axis=2))**2*\
            (np.sum(pi[np.newaxis,:,:]*pi[:,np.newaxis,:],axis=2))**2*\
            J_inter*np.exp(-rij_decay*(rijs_mag-rij_0))
        for i in range(0,n-1): #no nearest neighbor interactions
            wij[i,i+1] = 0
            wij[i+1,i] = 0
        H += wij
    #print(np.amax(wij))
    #max_pos = np.unravel_index(wij.argmax(),wij.shape)
    #print(max_pos[0])
    #print(r[max_pos[0],:])
    eigvals, eigvecs = LA.eigh(H) #eigenvals & vectors are listed from low to high eigenvalues
    return(eigvecs,eigvals,h_off,wij,H)

def plot_DOS(eigvals,res,mo0,start,stop):
    eigvals = eigvals.flatten()
    hist, bin_edges = np.histogram(eigvals,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt = np.array([np.average(eigvals),np.std(eigvals)])
    #popt,pcov = curve_fit(gaussian,bins,hist)
    fig, ax = plt.subplots()
    ax.hist(eigvals, bins=res, density=True)
    ax.plot(bins, gaussian(bins,popt[0],popt[1]), 'r-')
    ax.text(bins[int(bins.size/3*2)],gaussian(bins[int(bins.size/3*2)],popt[0],popt[1])*1.1+.02,r'$\mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1]))
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('DoS (per eV)')
    ax.set_title('DOS MO{} trjs{}-{}'.format(mo0,start,stop))
    plt.savefig('DOS_MO{}_trjs{}-{}.png'.format(mo0,start,stop))
    out = ('DOSave = {:.4f}\nDOSstd = {:.4f}'.format(popt[0],popt[1]))
    return(out)

def plot_wij(wij,res,start,stop,neighbors):
    trjs = wij.shape[0]
    n = wij.shape[1]
    wij = np.sort(wij, axis=-1)
    wij = wij[:,:,(n-(neighbors-1)):n]
    wij = wij.flatten()*1000
    wij = wij[ (wij > 1)]
    #popt1 = np.array([np.average(wij),np.std(wij),(np.size(wij)/trjs)])
    hist, bin_edges = np.histogram(wij,bins=res,range=(0,100))
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    i = 0
    while (hist[i] > 0):
        i += 1
        if i == res:
            break
    hist = hist/trjs
    #popt,pcov = curve_fit(gaussian,bins[1:i],(hist[1:i]))
    fig, ax = plt.subplots()
    ax.plot(bins, hist, 'bs')
    #ax.plot(bins[1:i], gaussian(bins[1:i],popt[0],popt[1]), 'r-')
    #ax.text(bins[int(i/2)],0.5,r'$\mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1]))
    ax.set_xlabel('Through Space Coupling (meV)')
    ax.set_ylabel('Counts (per trj)')
    ax.set_title('Top {} Wijs >1 meV trjs{}-{}'.format(neighbors,start,stop))
    plt.savefig('Wijs{}_trjs{}-{}.png'.format(neighbors,start,stop))
    #print('Wij top {} (no 0s) ave: {}'.format(neighbors, popt1[0]))
    #print('Wij top {} (no 0s) stdev: {}'.format(neighbors, popt1[1]))
    #print('Wij top {} (no 0s) number: {}'.format(neighbors, popt1[2]))
    return

def gaussian(x,mu,sig):
    return (1/(sig*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu)/sig)**2)

def gaussian0(x,sig): #gaussian centered at 0, with area = 2
    return (2/(sig*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-0)/sig)**2)

def gaussian2(x,mu1,sig1,mu2,sig2,ratio):
    return (ratio/(sig1*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu1)/sig1)**2+(1-ratio)/(sig2*math.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu2)/sig2)**2)

def exponential(x,a,b):
    return (a*np.exp(-b*x))

def linear(x,a,b):
    return(a*x+b)

def plot_deltabbpot(bbpots,res,start,stop):
    n = bbpots.shape[1]
    trjs = bbpots.shape[0]
    delta = bbpots[:,1:]-bbpots[:,:(n-1)]
    delta = delta.flatten()
    delta = np.absolute(delta)
    delta_max = np.amax(delta)/(1)
    hist, bin_edges = np.histogram(delta,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt,pcov = curve_fit(gaussian0,bins,hist)
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(8,4)) 
    ax1.plot(bins, hist/hist[0], 'bs')
    ax1.plot(bins, gaussian0(bins,popt[0])/hist[0], 'r-')
    ax1.text(bins[int(bins.size/2)],gaussian0(bins[0],popt[0])/2/hist[0],r'$\sigma={:.3f} eV$'.format(popt[0]))
    ax1.set_xlabel('Delta Electrostatic Potential (eV)')
    ax1.set_ylabel('Histogram (arbitrary units)')
    ax1.set_title('Delta Electrostatic Potential\nAlong the Backbone')
    out1 = ('DeltaEPalongBBstd = {:.4f}\n'.format(popt[0]))
    ijs = int((n**2-n)/2) #number of unique ij combinations
    delta_ijs = np.zeros((trjs,ijs),dtype = float)
    ij = 0
    for i in range((n-1)): ###n-1
        delta_ij = np.zeros((trjs,(n-1-i)),dtype = float)
        delta_ij = bbpots[:,i,np.newaxis]-bbpots[:,(i+1):] #trj x n-i-1 matrix
        ijend = ij + (n-1-i)
        delta_ijs[:,ij:ijend] = delta_ij
        ij = ijend
    delta_ijs = delta_ijs.flatten()
    delta_ijs = np.absolute(delta_ijs)
    hist, bin_edges = np.histogram(delta_ijs,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt,pcov = curve_fit(gaussian0,bins,hist)
    ax2.plot(bins, hist/hist[0], 'bs')
    ax2.plot(bins, gaussian0(bins,popt[0])/hist[0], 'r-')
    ax2.text(bins[int(bins.size/2)],gaussian0(bins[0],popt[0])/2/hist[0],r'$\sigma={:.3f} eV$'.format(popt[0]))
    ax2.set_xlabel('Delta Electrostatic Potential (eV)')
    ax2.set_title('Delta Electrostatic Potential\nAcross the Backbone')
    plt.savefig('Delta_Electrostatic_Potential_trjs{}-{}.png'.format(start,stop))
    out2 = 'DeltaEPacrossBBstd = {:.4f}'.format(popt[0])
    #plt.show()
    out = out1+out2
    return(out)

def bbpot_histogram(bbpots,res,start,stop):
    n = bbpots.shape[1]
    trjs = bbpots.shape[0]
    aves = np.zeros((trjs),dtype = float)
    aves = np.average(bbpots,axis=1) #average ep for each trj
    bbpots = bbpots - aves[:,np.newaxis]
    bbpots = bbpots.flatten()
    hist, bin_edges = np.histogram(bbpots,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt = np.array([np.average(bbpots),np.std(bbpots)])
    #popt,pcov = curve_fit(gaussian,bins,hist)
    fig, ax = plt.subplots()
    ax.plot(bins, hist, 'bs')
    ax.plot(bins, gaussian(bins,popt[0],popt[1]), 'r-')
    ax.text(bins[int(bins.size/3*2)],gaussian(bins[int(bins.size/3*2)],popt[0],popt[1])*1.1+.02,r'$\mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1]))
    ax.set_xlabel('Potential (eV) normalized so average = 0')
    ax.set_ylabel('Histogram')
    ax.set_title('BBPot Histogram trjs{}-{}'.format(start,stop))
    plt.savefig('BBPot_trjs{}-{}.png'.format(start,stop))
    out = 'EPave = {:.4f}\nEPstd = {:.4f}'.format(popt[0],popt[1])
    return(out)

def bbpot_histogram2(bbpots,res,start,stop):
    n = bbpots.shape[1]
    trjs = bbpots.shape[0]
    aves = np.zeros((trjs),dtype = float)
    aves = np.average(bbpots,axis=1) #average ep for each trj
    bbpots = bbpots - aves[:,np.newaxis]
    bbpots = bbpots.flatten()
    hist, bin_edges = np.histogram(bbpots,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt,pcov = curve_fit(gaussian2,bins,hist)
    #popt = np.array([np.average(bbpots),np.std(bbpots)])
    fig, ax = plt.subplots()
    ax.plot(bins, hist, 'bs')
    ax.plot(bins, gaussian2(bins,popt[0],popt[1],popt[2],popt[3],popt[4]), 'r-')
    ax.text(bins[0],gaussian2(bins[int(res/2)],popt[0],popt[1],popt[2],popt[3],popt[4])*1.1+.02,r'$\mu={:.2f},\ \sigma={:.2f}. \mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1],popt[2],popt[3]))
    ax.set_xlabel('Potential (eV) normalized so average = 0')
    ax.set_ylabel('Histogram')
    ax.set_title('BBPot Histogram double gaussiantrjs{}-{}'.format(start,stop))
    plt.savefig('BBPot_trjs{}-{}_2gaus.png'.format(start,stop))
    return

def print_h_off_stats(h_off,start,stop):
    out = ('Dihedralave = {:.4f}\nDihedralstd = {:.4f}'.format(np.average(h_off),np.std(h_off)))
    return(out)

def plot_IPR(eigvecs, eigvals,res,start,stop):
    n = eigvals.shape[1]
    trjs = eigvals.shape[0]
    IPRs = np.zeros((trjs*n,2), dtype = float)
    IPRs[:,0] = eigvals.flatten()
    eigvecs = eigvecs**4
    IPRs_temp = (np.sum(eigvecs,1))**(-1)
    IPRs[:,1] = IPRs_temp.flatten()
    ##Plot IPR vs Energy Heat Map
    fig, ax = plt.subplots()
    ax.scatter(IPRs[:,0],IPRs[:,1], marker="+", alpha=(0.01)) #100k = shading
    ax.set_ylim(bottom=0, top=60)
    ax.axvline(0, color="grey", linewidth = 0.5)
    ax.set_xticks(np.arange(21)-10)
    center = (np.amin(IPRs[:,0])+np.amax(IPRs[:,0]))/2
    ax.set_xlim((center-3.2,center+3.2))
    ax.set_yticks(np.arange(7)*10)
    ax.set_xlabel('State Energy (eV)')
    ax.set_ylabel('IPR')
    ax.set_title('IPR vs. Energy')
    plt.savefig('IPR_vs_Energy_trjs{}-{}.png'.format(start,stop))
    #plt.show()
    ##Plot Ave IPR vs Energy
    Emin = np.amin(IPRs[:,0])
    Emax = np.amax(IPRs[:,0])
    width = (Emax-Emin)/res
    Ebins = np.linspace(Emin,Emax,res+1)
    Ebins = Ebins[:res]+width/2
    Ecounter = np.zeros(res)
    Esum = np.zeros(res)
    IPRs[:,0] = np.floor((IPRs[:,0]-Emin)/width)
    for i in range (n*trjs):
        Ebin = int(IPRs[i,0])
        if Ebin == res:
            Ebin -= 1
        Ecounter[Ebin] += 1
        Esum[Ebin] += IPRs[i,1]
    for i in range (0,res):
        if Ecounter[i] == 0:
            Ecounter[i] = 1
    IPRave = Esum/Ecounter
    fig, ax = plt.subplots()
    ax.plot(Ebins,IPRave)
    ax.set_ylim(bottom=0, top=50)
    ax.axvline(0, color="grey", linewidth = 0.5)
    ax.set_yticks(np.arange(11)*5)
    ax.set_xlabel('State Energy (eV)')
    ax.set_ylabel('Average IPR (number of monomers)')
    ax.set_title('Average IPR vs State Energy')
    plt.savefig('AveIPR_vs_Energy_trjs{}-{}.png'.format(start,stop))
    #plt.show()
    return

def write_eigvals(eigvals, eigvals_name):
    f = open(eigvals_name, 'w')
    for i in range (eigvals.shape[0]):
        f.write('Step {}\nState\tEnergy(eV)\n'.format(i))
        for j in range (0,eigvals.shape[1]):
            f.write('{}\t'.format(j))
            f.write('{:.2f}\n'.format(eigvals[i,j]))
    f.close
    return

def write_eigvecs(eigvecs, eigvecs_name):
    f = open(eigvecs_name, 'w')
    for i in range (eigvecs.shape[0]):
        f.write('Step {}\nThese are the coefficients of each eigenvector at each site. Each column is an eigenvector.\nState1\tState2\tState3\tetc\n'.format(i))
        for j in range (0,eigvecs.shape[1]):
            for k in range (0,eigvecs.shape[2]): #this selects the vector
                f.write('{:.3f}\t'.format(eigvecs[i,j,k]))
            f.write('\n')
    f.close
    return

def write_bbpot(bbpots,bbpot_name):
    f = open(bbpot_name, 'w')
    for i in range (bbpots.shape[0]):
        f.write('Step {}\nSite\tElectrostatic Potential(eV)\n'.format(i))
        for j in range (0,bbpots.shape[1]):
            f.write('{}\t'.format(j))
            f.write('{:.4f}\n'.format(bbpots[i,j]))
    f.close
    return

def write_dihedrals(bbpots,bbpot_name):
    f = open(bbpot_name, 'w')
    for i in range (bbpots.shape[0]):
        f.write('Step {}\nSite\tDihedral Angles (degrees)\n'.format(i))
        for j in range (0,bbpots.shape[1]):
            f.write('{}\t'.format(j))
            f.write('{:.4f}\n'.format(bbpots[i,j]))
    f.close
    return

def plot_dihedrals(dihedrals,start,stop):
    boltzman = np.zeros((37,2),dtype = float)
    k_Angle = 2.4
    for i in range(37):
        angle = (i-18)*(math.pi/18)
        boltzman[i,0] = angle
        boltzman[i,1] = math.exp(-0.5*k_Angle*(1-math.cos(2*angle)))
    partition_function = np.sum(boltzman[:,1],axis=0)-boltzman[0,1] #don't double count
    boltzman[:,1] = boltzman[:,1]/partition_function*dihedrals.size
    fig, ax = plt.subplots()
    ax.plot(boltzman[:,0],boltzman[:,1],label='Boltzman')
    ax.hist(dihedrals.flatten(), bins=36, range=[-3.15,3.15], log=False, label='Histogram')
    ax.legend()
    ax.set_xlabel('Dihedral (radians)')
    ax.set_ylabel('Number of counts')
    ax.set_title('Histogram of Dihedrals compared to Boltzman Distribution')
    plt.savefig('Dihedrals_trjs{}-{}.png'.format(start,stop))
    #plt.show()
    return()

def plot_MOgaps(eigvals,res,start,stop,MOs):
    gaps = np.zeros((eigvals.shape[0],eigvals.shape[1]-1), dtype = float)
    gaps = eigvals[:,-1,np.newaxis]-eigvals[:,:-1]
    fig, ax = plt.subplots()
    for i in range(1,MOs):
        ax.hist(gaps[:,-i].flatten(), bins=res, density=True, histtype='step', label='HOMO-{}'.format(i))
    ax.legend()
    ax.set_xlim((-0.02,0.62))
    ax.set_xlabel('Energy Gap (eV)')
    ax.set_ylabel('Histogram (density)')
    ax.set_title('Histogram of MO Energy Gaps')
    plt.savefig('MOgaps{}-{}.png'.format(start,stop))
    #plt.show()
    return

def plot_MOenergies(eigvals,res,start,stop,MOs):
    fig, ax = plt.subplots()
    for i in range(1,MOs+1):
        ax.hist(eigvals[:,-i].flatten(), bins=res, density=True, histtype='step', label='HOMO-{}'.format(i-1))
    ax.legend()
    maxev = np.amax(eigvals)
    ax.set_xlim((maxev-0.2,maxev+.02))
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Histogram (density)')
    ax.set_title('Histogram of MO Energies')
    plt.savefig('MOenergies{}-{}.png'.format(start,stop))
    #plt.show()
    return

def plot_HOMOvsgap(eigvals,start,stop,MO):
    energies = eigvals[:,-1]
    gaps = eigvals[:,-1]-eigvals[:,(-1-MO)]
    fig, ax = plt.subplots()
    ax.scatter(energies,gaps, marker="+", alpha=(0.3))
    ax.set_ylim(bottom=-0.05)
    ax.set_xlabel('HOMO Energy (eV)')
    ax.set_ylabel('HOMO-{} Energy Gap (eV)'.format(MO))
    ax.set_title('HOMO vs HOMO-{} Energy Gap'.format(MO))
    plt.savefig('HOMOvsGap{}_trj{}-{}.png'.format(MO,start,stop))
    #plt.show()
    return

def plot_MOIPRs(eigvecs, eigvals,res,start,stop, MOs):
    n = eigvals.shape[1]
    trjs = eigvals.shape[0]
    IPRs = np.zeros((trjs,n), dtype = float)
    eigvecs = eigvecs**4
    IPRs = (np.sum(eigvecs,1))**(-1)
    IPRsx = IPRs[:,-MOs:].flatten()
    hist, bin_edges = np.histogram(IPRsx,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    #popt,pcov = curve_fit(gaussian,bins,hist)
    popt = np.array([np.average(IPRsx),np.std(IPRsx)])
    fig, ax = plt.subplots()
    ax.plot(bins, hist, 'bs')
    ax.plot(bins, gaussian(bins,popt[0],popt[1]), 'r-')
    ax.text(bins[int(bins.size/3)],gaussian(bins[int(bins.size/3)],popt[0],popt[1])*1.1,r'$\mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1]))
    ax.set_xlim(0,100)
    ax.set_xlabel('IPR (sites)')
    ax.set_ylabel('Histogram')
    ax.set_title('Histogram of top {} MOs IPRs'.format(MOs))
    plt.savefig('HOMO-(HOMO-{})IPRs{}-{}.png'.format(MOs,start,stop))
    out = 'HOMOIPRave = {:.4f}\nHOMOIPRstd = {:.4f}'.format(popt[0],popt[1])
    #plt.show()    
    return(out)

def plot_MOIPRs_inrange(eigvecs, eigvals, res,start,stop,MO_range):
    n = eigvals.shape[1]
    trjs = eigvals.shape[0]
    IPRs = np.zeros((trjs,n), dtype = float)
    eigvecs = eigvecs**4
    IPRs = (np.sum(eigvecs,1))**(-1)
    for i in range(2): #trjs):
        homo = eigvals[i,-1]
        IPRs[i,:] = IPRs[i,:] *((homo-eigvals[i,:]) < MO_range)
    IPRsx = IPRs.flatten()
    IPRsx = IPRsx[IPRsx != 0]
    hist, bin_edges = np.histogram(IPRsx,bins=res,density=True)
    bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2
    popt = np.array([np.average(IPRsx),np.std(IPRsx)])
    #popt,pcov = curve_fit(gaussian,bins,hist)
    fig, ax = plt.subplots()
    ax.plot(bins, hist, 'bs')
    ax.plot(bins, gaussian(bins,popt[0],popt[1]), 'r-')
    ax.text(bins[int(bins.size/3)],gaussian(bins[int(bins.size/3)],popt[0],popt[1])*1.1,r'$\mu={:.2f},\ \sigma={:.2f}$'.format(popt[0],popt[1]))
    ax.set_xlim(0,100)
    ax.set_xlabel('IPR (sites)')
    ax.set_ylabel('Histogram')
    ax.set_title('Histogram of MO IPRs in {}eV'.format(MO_range))
    plt.savefig('MOIPRs_in{}eV{}-{}.png'.format(MO_range,start,stop))
    out = 'MOIPR-0.4ave = {:.4f}\nMOIPR-0.4std = {:.4f}'.format(popt[0],popt[1])
    return(out)


def plot_MOpositions(eigvecs,mo0,start,stop):
    eigvecs = eigvecs**2
    trjs = eigvecs.shape[0]
    ave = (np.sum(eigvecs,axis=0))/trjs
    fig, ax = plt.subplots()
    ax.bar(np.arange(np.shape(eigvecs)[1], dtype = int), ave)
    ax.set_xlabel('Monomer Site')
    ax.set_ylabel('Ave Probability')
    ax.set_title('Ave Site Contribution to MO{} trjs{}-{}'.format(mo0,start,stop))
    plt.savefig('MO{}positions_trjs{}-{}.png'.format(mo0,start,stop))
    #plt.show()    
    return

def write_trjwithMOs(trj_filename, eigvecs, mo0, mof, bbpots):
    f_in = open(trj_filename, 'r')
    f_out = open('ovitoMOs{}-{}bbpots.trj'.format(mo0,mof),'w')
    trjs = eigvecs.shape[0]
    MOs = eigvecs.shape[2]
    n = eigvecs.shape[1]
    atomEP = [0, 0, 2, -2]
    for i in range(trjs):
        lines = ''
        for j in range(8):
            line = f_in.readline()
            lines += line
        f_out.write(lines)
        line = f_in.readline()
        splitline = line.split()
        splitline.append('bbpot')
        for j in range(MOs):
            splitline.append('MO{}'.format(mo0+j))
        line = ''
        for k in range(len(splitline)):
            line += '{} '.format(splitline[k])
        line += '\n'
        f_out.write(line)
        for j in range(n):  #Note: the 1st atom must be a bb bead
            line = f_in.readline()
            splitline = line.split()
            splitline.append(bbpots[i,j])
            for k in range(MOs):
                splitline.append((eigvecs[i,j,k])**2)
            line = ''
            for k in range(len(splitline)):
                line += '{} '.format(splitline[k])
            line += '\n'
            f_out.write(line)
            for k in range(4): #assumes there are 4 non-bb atoms between each bb
                line = f_in.readline()
                splitline = line.split() 
                splitline.append('{}'.format(atomEP[k]))
                for l in range(MOs):
                    splitline.append('0')
                line = ''
                for l in range(len(splitline)):
                    line += '{} '.format(splitline[l])
                line += '\n'
                f_out.write(line)
    f_in.close
    f_out.close
    return
"""
def writeout(out, outname):
    f = open(outname, 'w')
    for i in range(len(out)):
        f.write(str(out[i])+'\n')
        print(str(out[i]))
    f.close()
    return()

