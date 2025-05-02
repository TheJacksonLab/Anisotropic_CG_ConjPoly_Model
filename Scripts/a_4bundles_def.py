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


def get_Aij(chains,chainlength,r,box,A_cut): #Builds Adjacency matrix
    n = chains*chainlength
    A = np.zeros((n,n), dtype = float)
    rijs = np.zeros((n,n,3), dtype = float)
    rijs = r[np.newaxis,:,:]-r[:,np.newaxis,:]
    rijs = rijs / box
    rijs = rijs - np.rint(rijs) #minimal image #slow
    rijs = rijs * box
    rijs_mag = LA.norm(rijs,axis=-1)
    A = np.where(rijs_mag < A_cut, 1, 0) #adjacent = < 1.3 sigma
    return(A)

def get_bundles(A,exclude,maxblock,chainlength):
    n = A.shape[-1]
    bundles = np.zeros(n,dtype=float)
    for i in range(n):
        newnet = np.zeros(n,dtype=bool)
        newnet[i] = 1
        bundles[i], innet = findbundle(A,n,exclude,maxblock,chainlength,newnet)
    return(bundles)

def findbundle(A,n,exclude,maxblock,chainlength,newnet):
    avail = np.ones(n,dtype=bool)
    innet = newnet[:] + 0
    bundle = 0
    while np.sum(newnet) != 0:
        js = np.where(newnet == True)[0] #[new monoers]
        newnet = np.zeros(n,dtype=bool)
        for j in js:
            #Block out neighboring (+/-exclude)  monomers on the same chain
            for space in range(-exclude,exclude+1):
                if math.floor((j+space)/chainlength) == math.floor(j/chainlength) and (j+space)>=0 and (j+space)<n:
                    avail[j+space] = 0
        if len(js) > maxblock:
            js = enforcemaxblock(js, maxblock,chainlength)
        for j in js:
            newnet = (A[j] & avail) | newnet
        bundle += 1
        for j in range(len(js)-1):
            if (js[j+1]-js[j]) > 2:
                bundle += 1
        innet += newnet
    return(bundle,innet)

def enforcemaxblock(js, maxblock,chainlength): #Ensure there are no block of >3 consecutive monomers
    j = 0
    while j < (len(js)-maxblock):
        #is there a block of (maxblock+1) consecutive monomers?
        if (js[j+maxblock]-js[j]) == maxblock and math.floor(js[j+maxblock]/chainlength) == math.floor(js[j]/chainlength):
            #find the length of the current block
            k = j + 0
            block = np.array([],dtype=int)
            while k < (len(js)-1):
                if (js[k+1]-js[k]) == 1:
                    block = np.append(js[k],block)
                    k += 1
                else:
                    break
            block = np.append(js[k],block)
            block = np.flip(block)
            block_ave = math.floor(np.average(block))
            k = 0
            while k < len(block):
                if np.absolute(js[j+k] - block_ave) > 1:
                    js = np.delete(js,(j+k))
                    block = np.delete(block,k)
                    k -= 1
                k += 1
        j += 1
    return(js)

def get_intrainterRDFs(r,A,exclude,maxblock,chainlength,box,d_r,nbins):
    #r = [n,(xyz)], A = [n,n (Adjacency)]
    n = A.shape[-1]
    chains = int(n/chainlength)
    volume = box[0]*box[1]*box[2]
    rijs = np.zeros((n,n,3), dtype = float)
    rijs = r[np.newaxis,:,:]-r[:,np.newaxis,:]
    rijs = rijs / box
    rijs = rijs - np.rint(rijs) #minimal image #slow
    rijs = rijs * box #box is [3]
    rij_mag = np.sqrt(np.sum(rijs**2,axis = -1)) #sum xyz = [n,n (r)]  ##This step is slow
    histi = np.zeros((2,chains,int(nbins)), dtype = float) #[(intra,inter),chains,bins]
    for a in range(chains): #For each chain
        net = np.array([],dtype = int)
        for i in range(chainlength):
            newnet = np.zeros(n,dtype=bool)
            newnet[i+(chainlength*a)] = 1
            bundle, innet = findbundle(A,n,exclude,maxblock,chainlength,newnet)
            newmonomers = np.where(innet == True)[0] #[new monoers]
            net = np.append(net,newmonomers)
        net = np.unique(net) #list of monomers in bundled with chain 'a'
        n_net = int(len(net))
        notnet = np.arange(n,dtype = int)
        notnet = np.delete(notnet, net)
        n_notnet = int(len(notnet))
        intra_rij = np.delete(np.delete(rij_mag,notnet,axis=0),notnet,axis=1)
        inter_rij = np.delete(np.delete(rij_mag,notnet,axis=0),net,axis=1)
        histi[0,a], edges = np.histogram(intra_rij,bins=nbins,range=(0.0,(d_r*nbins))) #intra
        histi[0,a] = norm_hist(histi[0,a],volume,d_r,nbins,n_net,n)
        histi[1,a], edges = np.histogram(inter_rij,bins=nbins,range=(0.0,(d_r*nbins))) #inter
        histi[1,a] = norm_hist(histi[1,a],volume,d_r,nbins,n_net,n)
    RDFs = np.average(histi, axis = 1) #average over all chains
    return(RDFs)

def norm_hist(hist, volume, d_r, nbins, n1, n2): #normalizes the histogram
    hist = hist/(n1)/(4/3*np.pi)/(n2/volume)
    for i in range (0,nbins):
        hist[i] = hist[i]/(((i+1)*d_r)**3-(i*d_r)**3)
    return(hist)

def get_CorlLeng(all_hist): #intersection of two lines
    difference = all_hist[1]-all_hist[2] #g_intra - g_inter
    i = 0
    while difference[i] > 0:
        i += 1
    CorlLeng = (all_hist[1,i-1]-all_hist[2,i-1])/((all_hist[1,i-1]-all_hist[2,i-1])+(all_hist[2,i]-all_hist[1,i]))*(all_hist[0,i]-all_hist[0,i-1])+all_hist[0,i-1]
    return(CorlLeng)

def plot_RDFs(hist,xmax,start,stop,CorlLeng): #plots the pair dist. function
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"
    fig, ax = plt.subplots(figsize=(4,4))
    colors = np.array([(0.1,0.9,0.1),(0.4,0.7,0.4),(0,0,0)],dtype = float)
    ax.loglog(hist[0], hist[1], color=colors[0], label = '$g_{intra}$') 
    ax.loglog(hist[0], hist[2], color=colors[1], label = '$g_{inter}$') 
    ax.loglog(hist[0], hist[3], color=colors[2], label = '$g_{total}$') 
    ax.set_ylabel('RDF', fontsize='x-large')
    ax.set_xlabel('Distance ($\sigma$)', fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.axhline(y=1, color="gray")
    ax.text(3,2, 'Correlation Length = {:.2f} $\sigma$'.format(CorlLeng),zorder=10,fontsize=10)
    ax.set(xlim=(2,xmax), ylim=(0.1,20))
    #labels = ax.get_xticklabels()
    #plt.setp(labels, rotation=45)
    plt.xticks([])
    plt.xticks([2,3,4,6,10,20],['2','3','4','6','10','20'])
    #locator = ax.get_major_locator()
    #ax.xaxis.set_major_locator(locator)
    #ticker.FixedLocator([2, 5, 10, 20])
    ax.legend()
    fig.tight_layout()
    plt.savefig('CorlLeng_trj{}-{}.png'.format(start,stop), dpi=200)
    return

def write_trjwithbundles(trj_filename, bundles,fileout,start):
    f_in = open(trj_filename, 'r')
    f_out = open(fileout,'w')
    trjs = bundles.shape[0]
    n = bundles.shape[1]
    atomEP = [0, 0, 0, 0]
    for i in range(start*(9+n*5)):
        line = f_in.readline()
    for i in range(trjs):
        lines = ''
        for j in range(8):
            line = f_in.readline()
            lines += line
        f_out.write(lines)
        line = f_in.readline()
        splitline = line.split()
        splitline.append('bundle')
        line = ''
        for k in range(len(splitline)):
            line += '{} '.format(splitline[k])
        line += '\n'
        f_out.write(line)
        for j in range(n):  #Note: the 1st atom must be a bb bead
            line = f_in.readline()
            splitline = line.split()
            splitline.append(bundles[i,j])
            line = ''
            for k in range(len(splitline)):
                line += '{} '.format(splitline[k])
            line += '\n'
            f_out.write(line)
            for k in range(4):
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

def bundle_histogram(bundles, trj_int,start,stop):
    bin_edges = np.arange(1,11,1)
    rawhist,bins = np.histogram(bundles, bins=bin_edges)
    ave_raw,std_raw = np.average(bundles), np.std(bundles)
    weightedhist,bins = np.histogram(bundles, bins=bin_edges, weights = 1/bundles)
    #ave_wei =
    width = 0.3
    fig, ax = plt.subplots()
    ax.bar(bin_edges[:-1]-width/2,rawhist,width,label='Raw Histogram (per monomer)')
    ax.bar(bin_edges[:-1]+width/2,weightedhist,width,label='Weighted Histogram (per bundle)')
    ax.set_xlabel('Bundle Size')
    ax.set_ylabel('number of monomers')
    ax.text(3,rawhist[2],'Monomer Ave = {:.2f}\nMonomer Std = {:.2f}'.format(ave_raw,std_raw))
    ax.legend()
    plt.savefig('BundleHistogram{}-{}.png'.format(start,stop))
    return(ave_raw,std_raw)

def bundle_histogram_time(bundles, trj_int,start,stop):
    maxbundle = int(np.amax(bundles))
    bin_edges = np.arange(1,maxbundle,1)
    trjs = int(stop-start)
    rawhist = np.zeros((trjs,int(maxbundle-2)),dtype = int)
    for i in range(stop-start):
        rawhist[i],bins = np.histogram(bundles[i], bins=bin_edges)
    time = np.arange(trjs)/trjs
    fig, ax = plt.subplots()
    for i in range(maxbundle-2):
        ax.plot(time+i+0.5,rawhist[:,i])
    ax.set_xlabel('Bundle Size')
    ax.set_ylabel('number of monomers')
    plt.savefig('BundleHistogramTime{}-{}.png'.format(start,stop))
    return()


def bundle_time_corl(data,trj_int,corl_max):
    # data is a 3D matrix [trjs,n,(n parameters)] 
    corl_max += 1
    trjs = data.shape[0]
    monos = data.shape[1]
    n = 1 # n = number of parameters
    corl = np.zeros((corl_max,monos,2), dtype = float)
    average = np.full_like(data,np.average(data),dtype=float)
    data -= average
    corl[0,:,0] = 0
    corl[0,:,1] = np.average((data[:,:]**2),axis = 0)
    for i in range(1,corl_max):
        corl[i,:,0] = trj_int*i
        corl[i,:,1] = np.average((data[:(-i),:]*data[i:,:]),axis = 0)
    ave_time_corl = np.average(corl, axis = 1)
    ave_time_corl[:,1] = ave_time_corl[:,1]/ave_time_corl[0,1]
    return(ave_time_corl)

def plot_bundle_time_corl(ave_time_corl,start,stop,corl_range_i,corl_range_f):
    #plotting averages across all chains
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"
    ln_corl = np.zeros_like(ave_time_corl[(corl_range_i-1):(corl_range_f),:], dtype = float)
    ln_corl[:,:] = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    ln_corl[:,1] = np.log(ln_corl[:,1]) #take the natural log of the distance corellation
    corl_time_fit = np.polynomial.polynomial.polyfit(ln_corl[:,0],ln_corl[:,1],1)
    T_cor = -1/corl_time_fit[1]
    corl_range = ave_time_corl[(corl_range_i-1):(corl_range_f),:]
    fig, ax = plt.subplots()
    time = ave_time_corl[:,0]
    fig, ax = plt.subplots()
    ax.plot(time,ave_time_corl[:,1], label='Bundle Size')
    ax.plot(corl_range[:,0],corl_range[:,1],'b+', label='T_cor fit data')
    #ax.plot(time,np.exp(-time/T_cor),'g-', label='Correlation fit')
    ax.annotate('T_cor = {:.2f}'.format(T_cor),(corl_range[0,0],corl_range[0,1]*0.3))
    ax.set_xlabel('Time Interval (tau)', fontsize='x-large')
    ax.set_ylabel('Normalized Autocorrelation Function', fontsize='x-large')
    ax.legend()
    if np.amin(ave_time_corl[:,1]) > 0:
        ax.set_ylim(bottom=0)
    #plt.title('Average Bundle Size Time Correlation')
    ax.axhline(0, color="grey", linewidth = 0.5)
    plt.savefig('Bundle_Time_Corl{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close
    return(T_cor)

def writeout(out, outname):
    f = open(outname, 'w')
    for i in range(len(out)):
        f.write(str(out[i])+'\n')
        print(str(out[i]))
    f.close()
    return()
