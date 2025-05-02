import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import sys
from a_4lammps import  *

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

runs = 5 #number Lammps runs
time_step, trjs = get_timestep('savefile_run.out',runs)
energies = get_data('savefile_run.out',trjs,time_step,runs) #time is in units of tau

#############################################################
# Output Options: energies #
#############################################################

#1st region is plotted, all regions are averaged
starts = [0  ]
stops = [500]
out = []

for i in range (0,len(starts)):
    start = starts[i] #the 1st thermo point it will plot from 
    stop = stops[i] #the last thermo point it will plot 
    averaging = 10 #number of thermos averaged together
    if i == 0:    
        plot_energies(energies[start:stop,:],averaging,start,stop)
    out.append(print_ave_energy(energies[start:stop,:],start,stop))
writeout(out, 'zout_lammps.out')
