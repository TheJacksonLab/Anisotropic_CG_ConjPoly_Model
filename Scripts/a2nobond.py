import numpy as np
#import matplotlib.pyplot as plt
import math
import sys
import pathlib
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

#sys.argv[1] is the run name (i.e. P3ImHT3)
file_name = 'D10GB10_20wv.data'

def get_data(file1,atoms):  #gets distance and potential
    line = file1.readline()
    data = np.zeros((atoms,12), dtype = 'object')
    i = 0
    while(line != '\n'):
        data[i,:] = line.split()
        line = file1.readline()
        i += 1
    output = np.zeros((atoms,7), dtype = 'object')
    output[:,0] = data[:,0]
    output[:,1] = data[:,5]
    output[:,2] = data[:,1]
    output[:,3] = data[:,6]
    output[:,4] = data[:,2]
    output[:,5] = data[:,3]
    output[:,6] = data[:,4]
    return(output)

def get_header(file1):
    header = ''
    i = 0
    line = file1.readline()
    header = line + '\n'
    while('Atoms' not in line):
        if 'atoms' in line:
            x = line.split()
            atoms = int(x[0])
        if 'atom' in line:
            line = line + '\n'
            header += line
        if 'hi' in line:
            header += line
        if i > 100: #the header should be < 1000 lines into the file
            print("Header not found!")
            exit()
        i += 1
        line = file1.readline()
    header += '\n' + 'Atoms\n' + '\n'
    file1.readline()
    return(header, atoms)

def get_dielectric():
    file4 = open('aP3ImHT.initialize','r')
    line = file4.readline()
    while('dielectric' not in line):
        line = file4.readline()
    x = line.split()
    dielectric = x[1]
    return(dielectric)

file1 = open(file_name, 'r')
header, atoms = get_header(file1)
output = get_data(file1, atoms)
file1.close()
out_name = 'bNoBond.data'
file2 = open(out_name, 'w')
file2.write(header)
for i in range(0,atoms):
    file2.write('{} {} {} {} {} {} {}\n'.format(output[i,0],output[i,1],output[i,2],output[i,3],output[i,4],output[i,5],output[i,6]))
file2.close()

dielectric = get_dielectric()

out_name = 'bP3ImHT_compute.initialize'
file3 = open(out_name, 'w')
file3.write('#Header\n\
\n\
dimension 3\n\
units lj\n\
boundary p p p\n\
atom_style full\n\
pair_style coul/long 5.0\n\
kspace_style pppm 1.0e-4\n\
\n\
read_data "bNoBond.data"\n\
\n\
mass * 1\n\
pair_coeff * *\n\
dielectric {}\n\
\n\
group ell type 1\n\
group anions type 5\n\
group real type 1 3 4 5\n\
\n\
set group anions charge -1.000000001\n\
set group ell charge 0.000000001\n\
compute pea real pe/atom pair kspace\n\
dump 4 ell custom 1 nbbpots.out id c_pea\n\
dump_modify 4 sort id\n\
rerun ovito.trj dump x y z\n\
undump 4'.format(dielectric))
file3.close()
