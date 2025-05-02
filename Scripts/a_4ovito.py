import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import time
from scipy.interpolate import interpn
from numpy import linalg as LA
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

############ Process Data ###############
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

def get_timestep(file_name):
    file1 = open(file_name,'r')
    line = file1.readline()
    i = 0
    while( "  Time step     : " not in line):
        line = file1.readline()
        i = i+1
        if i > 1000: #the header should be < 1000 lines into the file
            print("In get_timestep: Time step not found!")
            return(0)
    line = line.split()
    timestep = float(line[3])
    return(timestep)

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

def get_anions(ovito):
    nreal = ovito.shape[1]
    trjs = ovito.shape[0]
    nbb = 0
    for j in range(nreal):
        if ovito[1,j,1] == 5: #if atom is an anion
            nbb += 1
    anions = np.zeros((trjs,nbb,12),dtype=float)
    for i in range(trjs):
        k = 0
        for j in range(nreal):
            if ovito[i,j,1] == 5: #if atom is an anion
                anions[i,k,:] = ovito[i,j,:12]
                k += 1
    return(anions)

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

""" Not used replaced with vectorized forms
def rot_mat(quat):
	b = quat[5]
	c = quat[6]
	d = quat[7]
	a = quat[8]
	res = np.zeros((3,3))
	res[0][0] = a*a + b*b - c*c - d*d
	res[0][1] = 2*b*c - 2*a*d
	res[0][2] = 2*b*d + 2*a*c
	res[1][0] = 2*b*c + 2*a*d
	res[1][1] = a*a - b*b + c*c - d*d
	res[1][2] = 2*c*d - 2*a*b
	res[2][0] = 2*b*d - 2*a*c
	res[2][1] = 2*c*d + 2*a*b
	res[2][2] = a*a - b*b - c*c + d*d
	return res
    
def new_dihedral(b0,b1,b2): #not used anymore; incorporated into calc_dihedrals
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y,x)
"""

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

########### Dihedrals #############
    
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

def write_dihedrals(bbpots,bbpot_name):
    f = open(bbpot_name, 'w')
    for i in range (bbpots.shape[0]):
        f.write('Step {}\nSite\tDihedral Angles (radians)\n'.format(i))
        for j in range (0,bbpots.shape[1]):
            f.write('{}\t'.format(j))
            f.write('{:.4f}\n'.format(bbpots[i,j]))
    f.close
    return

def plot_dihedrals(dihedrals,start,stop,temp):
    boltzman = np.zeros((37,2),dtype = float)
    k_Angle = 2.4
    for i in range(37):
        angle = (i-18)*(math.pi/18)
        boltzman[i,0] = angle
        boltzman[i,1] = math.exp(-0.5*k_Angle/temp*(1-math.cos(2*angle)))
    partition_function = np.sum(boltzman[:,1],axis=0)-boltzman[0,1] #don't double count
    boltzman[:,1] = boltzman[:,1]/partition_function*dihedrals.size
    fig, ax = plt.subplots()
    ax.plot(boltzman[:,0],boltzman[:,1],label='Boltzman')
    ax.hist(dihedrals.flatten(), bins=36, range=[-3.15,3.15], log=False, label='Histogram')
    ax.legend()
    ax.set_xlabel('Dihedral (radians)')
    ax.set_ylabel('Number of counts')
    ax.set_title('Histogram of Dihedrals compared to Boltzman Distribution')
    plt.savefig('Dihedrals{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return()

def calc_dihedral_time_corl(data,corl_max,trj_int):
    # data is a 2D matrix, Time x n Dihedrals
    corl_max += 1
    trjs = data.shape[0]
    n = data.shape[1] # n = number of bonds
    xy_vectors = np.zeros((trjs,n,2), dtype = float)
    xy_vectors[:,:,0] = np.cos(data[:,:])
    xy_vectors[:,:,1] = np.sin(data[:,:])
    corl = np.zeros((corl_max,2), dtype = float)
    corl[0,1] = 1
    for i in range(1,corl_max):
        dots = np.zeros((trjs-i,n), dtype = float)
        dots = np.sum((xy_vectors[:(-i),:,:]*xy_vectors[i:,:,:]),axis = 2)
        corl[i,1] = np.average(dots)
        corl[i,0] = trj_int*i
    return(corl)

def plot_dihedral_time_corl(time_corl,start,stop):
    time = time_corl[:,0]
    fig, ax = plt.subplots()
    ax.plot(time,time_corl[:,1], label='Dihedral Time Correlation')
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Dihedral Correlation')
    if np.amin(time_corl[:,1]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('Dihedral Time Correlation')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('Dihed_Time_Corl{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return()

###################### Bond Correlation (L_p) ###############

def plot_bondcorrelation(ellipsoids,trj_i,trj_f,corl_max,corl_fit_i,corl_fit_f):
    corl_range_i = corl_fit_i #first correlation length to use in calculating persistence length, inclusive
    corl_range_f = corl_fit_f #final correlation length to use in calculating pers length, inclusive
    trjs = ellipsoids.shape[0]
    chains = ellipsoids.shape[1]
    n = ellipsoids.shape[2]
    # r = np.zeros(((trj_f-trj_i),chains,n,3), dtype = float)
    r = ellipsoids[:,:,:,2:5] #r only includes relevant trjs and xyz data
    corl = np.zeros((corl_max+1,2), dtype = float)
    corl[:,0] = np.arange(0, int(corl_max+1), dtype=int)
    b = r[:,:,1:,:]-r[:,:,:(-1),:] #b is the bond vectors
    norm = LA.norm(b, axis=-1)
    b = b/np.repeat(norm[:,:,:,np.newaxis],3,axis=3) #b's vectors are normalized
    corl[0,1] = 1
    for i in range(1,corl_max+1):
        dots = np.zeros((trjs,chains,n-1-i), dtype = float)
        dots = np.sum(b[:,:,:(-i),:]*b[:,:,i:,:],axis = 3)
        corl[i,1] = np.average(dots)
    for i in range(corl_range_f): #cut off the fit before when the correlation goes negative 
        if corl[i,1] <= 0.01:
            corl_range_f = i
            break
    ln_corl = np.zeros_like(corl[(corl_range_i-1):(corl_range_f),:], dtype = float)
    ln_corl[:,:] = corl[(corl_range_i-1):(corl_range_f),:]
    ln_corl[:,1] = np.log(ln_corl[:,1]) #take the natural log of the distance corellation
    pers_length_fit = np.polynomial.polynomial.polyfit(ln_corl[:,0],ln_corl[:,1],1)
    Lp = -1/pers_length_fit[1]
    corl_range = corl[(corl_range_i-1):(corl_range_f),:]
    fig, ax = plt.subplots()
    ax.plot(corl[:,0],corl[:,1], 'ro', label='Correlation Data')
    ax.plot(corl_range[:,0],corl_range[:,1],'b+', label='Data used to calulate persistence length')
    ax.plot(corl[:,0],np.exp(-corl[:,0]/Lp),'g-', label='Correlation fit')
    ax.annotate('Lp = {:.2f}'.format(Lp),(corl_range[0,0],corl_range[0,1]))
    ax.legend()
    ax.set_title('Bond-Bond Correlation for trjs {}-{}'.format(trj_i,trj_f))
    ax.set_xlabel('Separation (BB bonds)')
    ax.set_ylabel('Bond Correlation')
    ax.set_ylim([-1,1])
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('BondCorrelation{}_{}.png'.format(trj_i,trj_f))
    #plt.show()
    plt.close
    return(Lp,int(corl_range_f-corl_range_i))

################ Order Parameters ################

def plot_axis_corl(ellipsoids,trj_i,trj_f,corl_max):
    trjs = ellipsoids.shape[0]
    n = ellipsoids.shape[1]
    r = np.zeros(((trj_f-trj_i),n,3), dtype = float)
    r = ellipsoids[trj_i:trj_f,:,2:5] #r only includes relevant trjs and xyz data
    corl = np.zeros((corl_max+1,2), dtype = float)
    corl[:,0] = np.arange(0, int(corl_max+1), dtype=int)
    b = r[:,1:,:]-r[:,:(-1),:] #b is the bond vectors
    norm = LA.norm(b, axis=-1)
    b = b/np.repeat(norm[:,:,np.newaxis],3,axis=2) #b's vectors are normalized
    #cross product of b1 and b2
    axes = np.zeros((trjs,(n-2),3),dtype = float)
    axes = np.cross(b[:,1:,:],b[:,:(-1),:])
    for i in range(0,1):
        dots = np.zeros((trjs,n-2-i), dtype = float)
        dots = np.sum(axes[:,:,:]*axes[:,:,:],axis = 2)
        #should I take the absolute value? kinks will reverse the axis
        corl[i,1] = np.average(dots)
    for i in range(1,corl_max+1):
        dots = np.zeros((trjs,n-2-i), dtype = float)
        dots = np.sum(axes[:,:(-i),:]*axes[:,i:,:],axis = 2)
        #should I take the absolute value? kinks will reverse the axis
        corl[i,1] = np.average(dots)
    corl[:,1] = corl[:,1]/corl[0,1] #norm all values by the 1st
    fig, ax = plt.subplots()
    ax.plot(corl[:,0],corl[:,1], 'ro', label='Correlation Data')
    ax.legend()
    ax.set_title('Axis-Axis Correlation for trjs {}-{}'.format(trj_i,trj_f))
    ax.set_xlabel('Separation (BB bonds)')
    ax.set_ylabel('Axis Correlation')
    ax.set_ylim([-1,1])
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('AxisCorrelation{}_{}.png'.format(trj_i,trj_f))
    #plt.show()
    plt.close
    return()

def plot_axis_op(ellipsoids,trj_i,trj_f,corl_max,d_r):
    trjs = trj_f-trj_i
    n = ellipsoids.shape[1]
    r = np.zeros(((trj_f-trj_i),n,3), dtype = float)
    r = ellipsoids[trj_i:trj_f,:,2:5] #r only includes relevant trjs and xyz data
    bins = math.floor(corl_max/d_r)
    hist = np.zeros(bins,dtype=float)
    histn = np.zeros(bins, dtype=int)
    b = r[:,1:,:]-r[:,:(-1),:] #b is the bond vectors
    norm = LA.norm(b, axis=-1)
    b = b/np.repeat(norm[:,:,np.newaxis],3,axis=2) #b's vectors are normalized
    #cross product of b1 and b2
    axes = np.zeros((trjs,(n-2),3),dtype = float)
    axes = np.cross(b[:,1:,:],b[:,:(-1),:]) #a trjx(n-2)x3 vector
    normarray = np.sum(axes[:,:,:]*axes[:,:,:],axis = 2) #all axis dotted with themselves
    normaxis = np.average(normarray)
    for i in range(trjs):
        dots = np.zeros((n-2,n-2), dtype = float)
        dots = np.sum(axes[i,:,np.newaxis,:]*axes[i,np.newaxis,:,:],axis = 2) #nxn matrix of dots
        dots = np.triu(dots,k=1) #a vector of the upper right triangle
        rij = np.zeros((n-2,n-2,3),dtype = float)
        rij = r[i,np.newaxis,1:-1,:]-r[i,1:-1,np.newaxis,:] # matrix of rijs
        rij_mag = np.sqrt(np.sum(rij**2,axis = -1)) #sum xyz
        rij_mag = np.triu(rij_mag,k=1) #a vector of the upper right triangle
        histi, binsi = np.histogram(rij_mag,bins=bins,range=(0,bins*d_r),weights=dots,density=False)
        histni, binsi = np.histogram(rij_mag,bins=bins,range=(0,bins*d_r),density=False)
        histi[0] = 0
        histni[0] = 0
        hist += histi
        histn += histni
        #if i%1000 == 0:
            #print("At {}/{} trjs".format(i,trjs))
    histn = np.where(histn < 1, 1 , histn)
    hist = hist/histn
    hist = hist/normaxis
    fig, ax = plt.subplots()
    ax.plot(binsi[1:],hist, 'ro', label='helix OP')
    ax.set_title('Helical OP trjs {}-{}'.format(trj_i,trj_f))
    ax.set_xlabel('Distance (sigma)')
    ax.set_ylabel('Helical OP')
    ax.axhline(0, color="grey", linewidth = 0.5)
    ax.set_ylim([-1.1,1.1])
    plt.savefig('HelicalOP{}_{}.png'.format(trj_i,trj_f))
    #plt.show()    
    plt.close
    return()

def plot_pinematic_op(vects,quats,start,stop,nematic_op_max,d_r):
    n = vects.shape[1]
    trjs = vects.shape[0]
    r = quats[:,:,2:5] #xyz of each ellipsoid
    bins = math.floor(nematic_op_max/d_r)
    hist = np.zeros(bins,dtype=float)
    histn = np.zeros(bins, dtype=int)
    for i in range(trjs): #this loop takes about 1 second per round because there are no neighbor lists used.
        rij = np.zeros((n,n,3),dtype = float)
        rij = r[i,np.newaxis,:,:]-r[i,:,np.newaxis,:] # n,n,3 matrix of rijs
        rij_mag = np.sum(np.square(rij),axis = -1) #sum xyz
        rij_mag = np.triu(rij_mag,k=1) #a vector of the upper right triangle
        rij_mag = np.sqrt(rij_mag) #n,n matrix
        #in_range = rij_mag < nematic_op_max
        #in_range = np.triu(in_range,k=1)
        dots = np.zeros((n,n))
        #dots = np.where(in_range, np.absolute(np.sum(vects[i,:,np.newaxis,:]*vects[i,np.newaxis,:,:],axis = 2)), 0) #nxn matrix of dots     #This line of code takes longer than just doing the sum for all points in the matrix... no idea why.
        dots = np.absolute(np.sum(vects[i,:,np.newaxis,:]*vects[i,np.newaxis,:,:],axis = 2)) #nxn matrix of dots
        dots = np.triu(dots,k=1) #a vector of the upper right triangle
        dots = 3/2*(dots**2)-1/2
        histi, binsi = np.histogram(rij_mag,bins=bins,range=(0,bins*d_r),weights=dots,density=False)
        histni, binsi = np.histogram(rij_mag,bins=bins,range=(0,bins*d_r),density=False)
        histi[0] = 0
        histni[0] = 0
        hist += histi
        histn += histni
        #if i%1000 == 0:
            #print("At {}/{} trjs".format(i,trjs))
    histn = np.where(histn < 1, 1 , histn)
    hist = hist/histn
    fig, ax = plt.subplots()
    ax.plot(binsi[1:],hist, 'ro', label='nematic OP')
    ax.legend()
    ax.set_title('Pi-Vector Nematic OP trjs {}-{}'.format(start,stop))
    ax.set_xlabel('Distance (sigma)')
    ax.set_ylabel('Nematic OP')
    ax.set_ylim([-0.5,1])
    plt.savefig('PiNematicOP{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return()
 
def nemOPs(pi,bond,start,stop):
    n = pi.shape[1]
    trjs = pi.shape[0]
    #axes = np.zeros((trjs,(n-2),3),dtype = float)
    #axes = np.cross(bond[:,1:,:],bond[:,:(-1),:]) #a trjx(n-2)x3 vector
    #norm = LA.norm(axes, axis=-1)
    #axes = axes/np.repeat(norm[:,:,np.newaxis],3,axis=2) #axis vectors are normalized
    #normarray = np.sum(axes[:,:,:]*axes[:,:,:],axis = 2) #all axis dotted with themselves
    #normaxis = np.average(normarray)
    Q_pi = calcQ(pi)
    Q_bond = calcQ(bond)
    #Q_axes = calcQ(axes)
    out = ('Q_pi = {:.4f}\nQ_bond = {:.4f}\n'.format(Q_pi,Q_bond))
    return(out)

def calcQ(vect): #vect is trjx'n'x3 (n can be anything)
    trjs = vect.shape[0]
    n = vect.shape[1]
    I = np.zeros((trjs,n,3,3),dtype = float)
    for i in range(3):
        I[:,:,i,i]= 1
    Q = 1/n*np.sum(3/2*vect[:,:,:,np.newaxis]*vect[:,:,np.newaxis,:]-1/2*I,axis=1)
    eigvals = np.zeros((trjs,3), dtype = float)
    eigvecs = np.zeros((trjs,3,3), dtype = float)
    for i in range (trjs):
        eigvals[i,:], eigvecs[i,:,:] = LA.eigh(Q[i,:,:])
    Q1 = np.average(eigvals[:,2]) # the 3rd eigenvalue is the largest (always)
    return(Q1)

###################### Radius of Gyration & Asphere ... ################

def get_Rgs(ellips): #[trjs,chains,n,12] #notation from https://en.wikipedia.org/wiki/Gyration_tensor
    r = ellips[:,:,:,2:5] #r only includes relevant trjs and xyz data
    #r = np.array([[[[1,0,0],[0,0,0],[0,1,0],[1,1,0],[0.5,0.5,0.5],[0.5,0.5,-0.5]]]])
    center = np.average(r[:,:,:,:],axis = -2)
    center = np.repeat(center[:,:,np.newaxis,:],ellips.shape[2], axis = -2)
    r = r - center
    rnrms = r[:,:,:,:,np.newaxis]*r[:,:,:,np.newaxis,:]
    Snm = np.average(rnrms, axis = 2) #average over atoms
    eigvals, eigvecs = LA.eigh(Snm) #eigvals are low to high [x**2, y**2, z**2]
    output = np.zeros((ellips.shape[0],ellips.shape[1],4),dtype = float)
    output[:,:,0] = (np.sum(eigvals,axis = -1))**(0.5) #Radius of Gyration
    output[:,:,1] = (1.5*eigvals[:,:,2]-0.5*output[:,:,0]**2) #Asphericity
    output[:,:,2] = eigvals[:,:,1]-eigvals[:,:,0] #Acylindricity
    output[:,:,3] = (output[:,:,1]**2+0.75*output[:,:,2]**2)/(output[:,:,0]**4) #Anisotropy
    Rgs_ave = np.average(output, axis = (0,1))
    Rgs_std = np.std(output, axis = (0,1))
    return(output, Rgs_ave, Rgs_std)

def calc_rg_time_corls(data,corl_max,trj_int):
    # data is a 3D matrix [trjs,chains,(n parameters)] 
    corl_max += 1
    trjs = data.shape[0]
    chains = data.shape[1]
    n = data.shape[2] # n = number of parameters
    corl = np.zeros((corl_max,chains,5), dtype = float)
    ###Update Averaging: Average over all chains###
    average = np.average(data,axis = (0,1))
    data -= np.repeat(average[np.newaxis,np.newaxis,:], trjs, axis = 0)
    corl[0,:,0] = 0
    corl[0,:,1:] = np.average((data[:,:,:]**2),axis = 0)
    for i in range(1,corl_max):
        corl[i,:,0] = trj_int*i
        corl[i,:,1:] = np.average((data[:(-i),:,:]*data[i:,:,:]),axis = 0)
    corl[:,:,1:] = corl[:,:,1:]/corl[0,:,1:]
    return(corl)

def plot_rg_time_corl(time_corl,start,stop,corl_range_i,corl_range_f):
    #plotting averages across all chains
    ave_time_corl = np.average(time_corl,axis = 1)
    ln_corl = np.zeros_like(ave_time_corl[(corl_range_i-1):(corl_range_f),:], dtype = float)
    ln_corl[:,:] = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    ln_corl[:,1] = np.log(ln_corl[:,1]) #take the natural log of the distance corellation
    corl_time_fit = np.polynomial.polynomial.polyfit(ln_corl[:,0],ln_corl[:,1],1)
    T_cor = -1/corl_time_fit[1]
    corl_range = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    fig, ax = plt.subplots()
    time = ave_time_corl[:,0]
    labels = ['Time','Rg', 'Asphericity', 'Acylindricity', 'Anisotropy']
    fig, ax = plt.subplots()
    for i in range(1,5):
        ax.plot(time,ave_time_corl[:,i], label=labels[i])
    ax.plot(corl_range[:,0],corl_range[:,1],'b+', label='Data used to calulate persistence length')
    ax.plot(time,np.exp(-time/T_cor),'g-', label='Correlation fit')
    ax.annotate('T_cor = {:.2f}'.format(T_cor),(corl_range[0,0],corl_range[0,1]))
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Normalized Autocorrelation Function')
    ax.legend()
    if np.amin(time_corl[:,1:5]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('Average Chain Shape Time Correlation')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('Rg_Time_Corl{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(T_cor)


################ End to End vector #######################

def get_ends_length_time_corl(ellipsoids,corl_max,trj_int):
    ends = ellipsoids[:,:,-1,2:5]-ellipsoids[:,:,0,2:5]
    norm = LA.norm(ends, axis = -1)
    corl_max += 1
    trjs = norm.shape[0]
    chains = norm.shape[1]
    average = np.average(norm)
    norm -= average
    corl = np.zeros((corl_max,chains,2), dtype = float)
    corl[0,:,0] = 0
    corl[0,:,1] = np.average((norm[:,:]**2),axis = 0)
    for i in range(1,corl_max):
        corl[i,:,0] = trj_int*i
        corl[i,:,1] = np.average((norm[:(-i),:]*norm[i:,:]),axis = 0)
    corl[:,:,1] = corl[:,:,1]/corl[0,:,1]
    return(corl) #[time_gap, chain, (time, e.e)]

def get_ends_vector_time_corl(ellipsoids,corl_max,trj_int):
    ends = ellipsoids[:,:,-1,2:5]-ellipsoids[:,:,0,2:5]
    norm = LA.norm(ends, axis = -1)
    ends = ends/np.repeat(norm[:,:,np.newaxis], 3, axis=-1) #normalized end-to-end vector
    corl_max += 1
    trjs = ends.shape[0]
    chains = ends.shape[1]
    corl = np.zeros((corl_max,chains,2), dtype = float)
    corl[0,:,0] = 0 #time
    corl[0,:,1] = 1 #corl
    for i in range(1,corl_max):
        corl[i,:,0] = trj_int*i
        dots = np.zeros((trjs-i,chains), dtype = float)
        dots = np.sum((ends[:(-i),:,:]*ends[i:,:,:]),axis = 2)
        corl[i,:,1] = np.average(dots,axis = 0)
    return(corl) #[time_gap, (time, 4 chain e.e)]

def plot_E2El_time_corl(time_corl,start,stop,corl_range_i,corl_range_f):
    #plotting averages across all chains
    ave_time_corl = np.average(time_corl,axis = 1)
    ln_corl = np.zeros_like(ave_time_corl[(corl_range_i-1):(corl_range_f),:], dtype = float)
    ln_corl[:,:] = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    ln_corl[:,1] = np.log(ln_corl[:,1]) #take the natural log of the distance corellation
    corl_time_fit = np.polynomial.polynomial.polyfit(ln_corl[:,0],ln_corl[:,1],1)
    T_cor = -1/corl_time_fit[1]
    corl_range = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    time = ave_time_corl[:,0]
    labels = ['Time','End-to-End Orientation']
    fig, ax = plt.subplots()
    for i in range(1,2):
        ax.plot(time,ave_time_corl[:,i], label=labels[i])
    ax.plot(corl_range[:,0],corl_range[:,1],'b+', label='Data used to calulate persistence length')
    ax.plot(time,np.exp(-time/T_cor),'g-', label='Correlation fit')
    ax.annotate('T_cor = {:.2f}'.format(T_cor),(corl_range[0,0],corl_range[0,1]))
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Autocorrelation Function')
    ax.legend()
    if np.amin(time_corl[:,1:5]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('End-to-End Length Time Correlation')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('E2El_Time_Corl{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(T_cor)

def plot_E2Ev_time_corl(time_corl,start,stop,corl_range_i,corl_range_f):
    #plotting averages across all chains
    ave_time_corl = np.average(time_corl,axis = 1)
    ln_corl = np.zeros_like(ave_time_corl[(corl_range_i-1):(corl_range_f),:], dtype = float)
    ln_corl[:,:] = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    ln_corl[:,1] = np.log(ln_corl[:,1]) #take the natural log of the distance corellation
    corl_time_fit = np.polynomial.polynomial.polyfit(ln_corl[:,0],ln_corl[:,1],1)
    T_cor = -1/corl_time_fit[1]
    corl_range = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    time = ave_time_corl[:,0]
    labels = ['Time','End-to-End Orientation']
    fig, ax = plt.subplots()
    for i in range(1,2):
        ax.plot(time,ave_time_corl[:,i], label=labels[i])
    ax.plot(corl_range[:,0],corl_range[:,1],'b+', label='Data used to calulate persistence length')
    ax.plot(time,np.exp(-time/T_cor),'g-', label='Correlation fit')
    ax.annotate('T_cor = {:.2f}'.format(T_cor),(corl_range[0,0],corl_range[0,1]))
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Autocorrelation Function')
    ax.legend()
    if np.amin(time_corl[:,1:5]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('End-to-End Vector Orientation Time Correlation')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('E2Ev_Time_Corl{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(T_cor)

################ Polymer Diffusion (Mean Squared Displacement) ###################

def calc_MeSqDisp(ellipsoids,MSD_max,trj_int):
    r = ellipsoids[:,:,:,2:5] #r only includes relevant trjs and xyz data
    COMs = np.average(r, axis = 2)
    MSD_max += 1
    trjs = COMs.shape[0]
    chains = COMs.shape[1]
    MSD = np.zeros((MSD_max,chains,5), dtype = float) #time, r, x, y, z
    MSD[0,:,0] = 0 #time
    MSD[0,:,1:] = 0 #MSD
    for i in range(1,MSD_max):
        MSD[i,:,0] = trj_int*i
        delta2 = np.zeros((trjs-i,chains), dtype = float)
        delta2 = (COMs[:(-i),:,:]-COMs[i:,:,:])**2
        MSD[i,:,1] = np.average(np.sum(delta2,axis=2),axis=0) #MSD_r
        MSD[i,:,2:] = np.average(delta2,axis=0) #MSD_x,y,z
    return(MSD)

def plot_MeSqDisp(MSD,start,stop,D_range_i,D_range_f):
    #plotting averages across all chains
    ave_MSD = np.average(MSD,axis = 1)
    fig, ax = plt.subplots()
    time = ave_MSD[:,0]
    labels = ['Time','Dist2', 'x2', 'y2', 'z2']
    D_range = ave_MSD[(D_range_i-1):(D_range_f),:]
    MSD_fit = np.polynomial.polynomial.polyfit(D_range[:,0],D_range[:,1],1)
    Diff = MSD_fit[1]/6 #In 3D, MSD = 6*D*time
    MSD_fit0 = np.polynomial.polynomial.polyfit(ave_MSD[:2,0],ave_MSD[:2,1],1)
    Diff0 = MSD_fit0[1]/6 #In 3D, MSD = 6*D*time
    fig, ax = plt.subplots()
    for i in range(1,5):
        ax.plot(time,ave_MSD[:,i], label=labels[i])
    ax.plot(D_range[:,0],D_range[:,1],'b+', label='Data used to calulate Diff')
    ax.plot(time,(6*time*Diff),'g-', label='Diffusion Fit')
    ax.annotate('Diff = {:.2f} sigma2/ktau'.format(Diff*1000),(D_range[0,0],D_range[0,1]))
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Mean Squared Displacement')
    ax.legend()
    if np.amin(ave_MSD[:,1:5]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('Average Mean Squared Displacement')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('Poly_MeSqDisp{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(Diff, Diff0)

################ Anion Diffusion (Mean Squared Displacement) ###################

def calc_Anion_MeSqDisp(anions,MSD_max,trj_int):
    r = anions[:,:,2:5] #r only includes relevant trjs and xyz data
    MSD_max += 1
    trjs = r.shape[0]
    anions = r.shape[1]
    MSD = np.zeros((MSD_max,anions,5), dtype = float) #time, r, x, y, z
    MSD[0,:,0] = 0 #time
    MSD[0,:,1:] = 0 #MSD
    for i in range(1,MSD_max):
        MSD[i,:,0] = trj_int*i
        delta2 = np.zeros((trjs-i,anions), dtype = float)
        delta2 = (r[:(-i),:,:]-r[i:,:,:])**2
        MSD[i,:,1] = np.average(np.sum(delta2,axis=2),axis=0) #MSD_r
        MSD[i,:,2:] = np.average(delta2,axis=0) #MSD_x,y,z
    return(MSD)

def plot_Anion_MeSqDisp(MSD,start,stop,D_range_i,D_range_f):
    #plotting averages across all chains
    ave_MSD = np.average(MSD,axis = 1)
    fig, ax = plt.subplots()
    time = ave_MSD[:,0]
    labels = ['Time','Dist2', 'x2', 'y2', 'z2']
    D_range = ave_MSD[(D_range_i-1):(D_range_f),:]
    MSD_fit = np.polynomial.polynomial.polyfit(D_range[:,0],D_range[:,1],1)
    Diff = MSD_fit[1]/6 #In 3D, MSD = 6*D*time
    MSD_fit0 = np.polynomial.polynomial.polyfit(ave_MSD[:2,0],ave_MSD[:2,1],1)
    Diff0 = MSD_fit0[1]/6 #In 3D, MSD = 6*D*time
    fig, ax = plt.subplots()
    for i in range(1,5):
        ax.plot(time,ave_MSD[:,i], label=labels[i])
    ax.plot(D_range[:,0],D_range[:,1],'b+', label='Data used to calulate Diff')
    ax.plot(time,(6*time*Diff),'g-', label='Anion Diffusion Fit')
    ax.annotate('Diff = {:.2f} sigma2/ktau'.format(Diff*1000),(D_range[0,0],D_range[0,1]))
    ax.set_xlabel('Time Interval (Tau)')
    ax.set_ylabel('Mean Squared Displacement')
    ax.legend()
    if np.amin(ave_MSD[:,1:5]) > 0:
        ax.set_ylim(bottom=0)
    plt.title('Average Mean Squared Displacement')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('Anion_MeSqDisp{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(Diff,Diff0)

##################### Write Out ###############
def writeout(out, outname):
    f = open(outname, 'w')
    for i in range(len(out)):
        f.write(out[i]+'\n')
        #print(out[i])
    f.close()
    return()
