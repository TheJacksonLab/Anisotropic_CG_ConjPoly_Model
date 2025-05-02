import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.interpolate import interpn
from numpy import linalg as LA
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

def get_timestep(file_name,runs):
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
    g = 0
    thermos = 0
    for g in range(0,runs):
        while(line != "Step Temp c_ellTemp Press Volume c_rg c_gbpair c_ljpair c_coulpair E_long E_pair E_bond E_angle PotEng KinEng TotEng \n"):
            line = file1.readline()
            i = i+1
            if i > 5000: #the header should be < 1000 lines into the file
                print("In get_timestep: Header not found! g={}".format(g))
                return()
        while ("Loop" not in line):
            thermos += 1
            line = file1.readline()
        g += 1
        i = 0
    return(timestep,thermos)

def get_data(file_name,trjs,time_step,runs):
    rgs = np.zeros((trjs,2), dtype = 'float')
    gyration = np.zeros((trjs,7), dtype = 'float')
    energies = np.zeros((trjs,11), dtype = 'float')
    file1 = open(file_name,'r')
    line = file1.readline()
    g = 0
    h = 0
    i = 0
    for g in range(0,runs):
        while(line != "Step Temp c_ellTemp Press Volume c_rg c_gbpair c_ljpair c_coulpair E_long E_pair E_bond E_angle PotEng KinEng TotEng \n"):
            line = file1.readline()
            h += 1
            if h > 1000: #the header should be < 1000 lines into the file
                print("In get_rgs: Header not found! h = {}".format(h))
                return()
        line = file1.readline()
        h = 0
        while "Loop" not in line:
            line = line.split()
            energies[i,0] = float(line[0])*time_step
            energies[i,1] = float(line[6])
            energies[i,2] = float(line[7])
            energies[i,3] = float(line[8])
            energies[i,4] = float(line[9])
            energies[i,5] = float(line[10])
            energies[i,6] = float(line[11])
            energies[i,7] = float(line[12])
            energies[i,8] = float(line[13])
            energies[i,9] = float(line[14])
            energies[i,10] = float(line[15])
            line = file1.readline()
            i += 1
            h += 1
            if h > 1000000: # there should not be > 1M thermos
                print("In get_rgs: End of thermos not found!!")
                return(energies)
        g += 1
        h = 0
    return(energies)

def plot_energies(energies,ave_interval,start,stop):
    ave_points = math.floor(energies.shape[0]/ave_interval)
    num_energies = energies.shape[1]
    ave_energies = np.zeros((ave_points,num_energies), dtype = float)
    for i in range(ave_points):
        ave_energies[i,:] = np.average(energies[ave_interval*i:(ave_interval*(i+1)),:],axis=0)
    time = ave_energies[:,0]/1000 #1000 converts from Tau to kTau
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.subplots_adjust(wspace=0.5)    
    ax2.plot(time,ave_energies[:,1], label='GB')
    ax2.plot(time,ave_energies[:,2], label='LJ')
    ax2.plot(time,ave_energies[:,3], label='Coul/Cut')
    ax2.plot(time,ave_energies[:,4], label='Coul/Long')
    ax1.plot(time,ave_energies[:,5], label='Pair (total)')
    ax1.plot(time,ave_energies[:,6], label='Bond')
    ax1.plot(time,ave_energies[:,7], label='Angle')
    ax1.plot(time,ave_energies[:,8], label='PE (total)')
    ax1.plot(time,ave_energies[:,9], label='KE')
    #ax1.plot(time,ave_energies[:,10], label='Total Energy')
    ax1.set_xlabel('Time (kTau)')
    ax1.set_ylabel('Energy (kT)')
    ax1.set_title('Energies vs Time')
    ax1.axhline(0, color='grey', linewidth=0.5)
    ax1.legend()
    ax2.set_xlabel('Time (kTau)')
    ax2.set_ylabel('Energy (kT)')
    ax2.set_title('Pair Energies')
    ax2.axhline(0, color='grey', linewidth=0.5)
    ax2.legend()
    plt.savefig('Energies{}_{}.png'.format(start,stop))
    #plt.show()
    plt.close()
    return()

def print_ave_energy(energies,start,stop):
    ave = np.average(energies,axis=0)
    std = np.std(energies,axis=0)
    out = ('Thermos {}-{}\n\
GBave = {:.4f}\n\GBstd = {:.4f}\n\
LJave = {:.4f}\n\LJstd = {:.4f}\n\
Coul/Cutave = {:.4f}\n\Coul/Cutstd = {:.4f}\n\
Coul/Longave = {:.4f}\n\Coul/Longstd = {:.4f}\n\
Pair(total)ave = {:.4f}\n\Pair(total)std = {:.4f}\n\
Bondave = {:.4f}\n\Bondstd = {:.4f}\n\
Angleave = {:.4f}\n\Anglestd = {:.4f}\n\
PE(total)ave = {:.4f}\n\PE(total)std = {:.4f}\n\
KEave = {:.4f}\n\KAstd = {:.4f}\n\
TotalEnergyave = {:.4f}\n\TotalEnergystd = {:.4f}\n'.format(start,stop,ave[1],std[1],ave[2],std[2],ave[3],std[3],ave[4],std[4],ave[5],std[5],ave[6],std[6],ave[7],std[7],ave[8],std[8],ave[9],std[9],ave[10],std[10]))
    return(out)


def writeout(out, outname):
    f = open(outname, 'w')
    for i in range(len(out)):
        print(out[i])
        f.write(out[i]+'\n')
    f.close()
    return()
