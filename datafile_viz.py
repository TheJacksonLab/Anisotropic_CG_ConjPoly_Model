#Script for generating gay berne chains 
#have angles (between monomers and bond sites) and the special dihedrals that Alec wrote

import math
import numpy as np
from scipy import stats

#VARIABLES THAT USER CAN SET
NUM_MOLECULES = 1
NUM_MONOMERS = 32
monomerOffset = 2.4
moleculeOffset = 15

#file to write the data to
dataFile = open("viz.data","w")

#Box dimensions
box_x = 100.0
box_y = 100.0
box_z = 100.0

# Read in data file and store in array
monomer_data_atoms = []
monomer_data_angles = []
monomer_data_bonds = []
monomer_data_rigid = []
monomer_disks = 0
connecting_bonds = 0
connecting_angles = 0

mon_data = open("v5_nosc.mon", "r")
line = mon_data.readline()
while(line != "Particles:\n"):
	line = mon_data.readline()
line = mon_data.readline()
while(line != "Angles:\n"):
	splitline = line.split()
	if(int(splitline[6]) == 1):
		monomer_disks += 1
	temp = [int(splitline[0]), int(splitline[1]), float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), int(splitline[6]), float(splitline[7]), float(splitline[8]), float(splitline[9]), float(splitline[10])]
	monomer_data_atoms.append(temp)
	line = mon_data.readline()
line = mon_data.readline()
while(line != "Dihedrals/Bonds:\n"):
	splitline = line.split()
	if(int(splitline[0]) == -1):
		connecting_angles += 1
	temp = [int(splitline[0]), int(splitline[1]), int(splitline[2]), int(splitline[3]), int(splitline[4])]
	monomer_data_angles.append(temp)
	line = mon_data.readline()
line = mon_data.readline()
while(line != "Rigid:\n"):
	splitline = line.split()
	if(int(splitline[0]) == -1):
		connecting_bonds += 1
	temp = [int(splitline[0]), int(splitline[1]), int(splitline[2]), int(splitline[3])]
	monomer_data_bonds.append(temp)
	line = mon_data.readline()
line = mon_data.readline()
while(line != ""):
	splitline = line.split()
	temp = [int(splitline[0]), int(splitline[1])]
	monomer_data_rigid.append(temp)
	line = mon_data.readline()

#Number of ellipsoids in whole simulation
numEllipsoids = NUM_MOLECULES*NUM_MONOMERS*monomer_disks

#Number of monomers in the whole simulation
numMonomers = NUM_MOLECULES*NUM_MONOMERS

#Number of atoms per monomer
numAtomsPerMonomer = len(monomer_data_atoms)

#Number of atoms per chain
numAtomPerChain = NUM_MONOMERS*len(monomer_data_atoms)

#Number of atoms in the whole simulation
numAtom = numAtomPerChain * NUM_MOLECULES

#number of bonds per monomer
numBondsPerMonomer = len(monomer_data_bonds)

#number of bonds per molecule
numBondsPerMolecule = NUM_MONOMERS*len(monomer_data_bonds) - connecting_bonds

#number of bonds in the whole simulation
numBonds = numBondsPerMolecule*NUM_MOLECULES

#number of angles per monomer
numAnglesPerMonomer = len(monomer_data_angles)

#number of angles per molecule
numAnglesPerMolecule = NUM_MONOMERS*len(monomer_data_angles) - connecting_angles

#number of angles in the whole simulation
numAngles = numAnglesPerMolecule*NUM_MOLECULES

#number of atom, bond, and angle types in the simulation
atomTypes = int(np.amax(np.array(monomer_data_atoms)[:,1]))
bondTypes = int(np.amax(np.array(monomer_data_bonds)[:,1]))
angleTypes = int(np.amax(np.array(monomer_data_angles)[:,1]))

#Writing header of datafile
dataFile.write("Lammps Description\n\n")
dataFile.write("\t" + str(numAtom) + "\tatoms\n")
dataFile.write("\t" + str(numBonds) + "\tbonds\n")
dataFile.write("\t" + str(numAngles) + "\tangles\n")
dataFile.write("\t" + str(numEllipsoids) + "\tellipsoids\n")
dataFile.write("\n")

dataFile.write("\t" + str(atomTypes) + "\tatom types\n")
dataFile.write("\t" + str(bondTypes) + "\tbond types\n")
dataFile.write("\t" + str(angleTypes) + "\tangle types\n")
dataFile.write("\n")

dataFile.write(str(-1.0*box_x/2) + " " + str(box_x/2) + " " + "xlo xhi\n")
dataFile.write(str(-1.0*box_y/2) + " " + str(box_y/2) + " " + "ylo yhi\n")
dataFile.write(str(-1.0*box_z/2) + " " + str(box_z/2) + " " + "zlo zhi\n")
dataFile.write("\n")

dataFile.write("Atoms\n\n")

for mol in range(NUM_MOLECULES):
	startX = 0.0
	startY = 0.0 + moleculeOffset*mol
	startZ = 0.0
	for i in range(NUM_MONOMERS):
		x = startX + monomerOffset*i
		y = startY
		z = startZ
		#ID type x y z moleculeID charge ellipsoidFLAG density
		for d in monomer_data_atoms:
			dataFile.write(str(d[0] + mol*numAtomPerChain + i*numAtomsPerMonomer) + " " + str(d[1]) + " " + str(d[6]) + " " + str(d[7]) + " " + str(d[2]+x) + " " + str(d[3]+y) + " " + str(d[4]+z) + "\n")


dataFile.write("\n")

dataFile.write("Bonds\n\n")

for i in range(NUM_MOLECULES):
	for j in range(NUM_MONOMERS):
		n = 0
		tmp = []
		start = i*numAtomPerChain + j*numAtomsPerMonomer
		for d in monomer_data_bonds:
			if(d[0] > -1):
				dataFile.write(str(d[0] + j*numBondsPerMonomer + i*numBondsPerMolecule) + " " + str(d[1]) + " " + str(d[2] + start) + " " + str(d[3] + start) + "\n")
				n = d[0] + j*numBondsPerMonomer + i*numBondsPerMolecule
			if(d[0] == -1):
				tmp.append(d)
		if(j != NUM_MONOMERS-1):
			for t in tmp:
				dataFile.write(str(n+1) + " " + str(t[1]) + " " + str(t[2] + start) + " " + str(t[3] + start) + "\n")
				n += 1
dataFile.write("\n")

dataFile.write("Angles\n\n")
for i in range(NUM_MOLECULES):
	for j in range(NUM_MONOMERS):
		n = 0
		tmp = []
		start = i*numAtomPerChain + j*numAtomsPerMonomer
		for d in monomer_data_angles:
			if(d[0] > -1):
				dataFile.write(str(d[0] + j*numAnglesPerMonomer + i*numAnglesPerMolecule) + " " + str(d[1]) + " " + str(d[2] + start) + " " + str(d[3] + start) + " " + str(d[4] + start) + "\n")
				n = d[0] + j*numAnglesPerMonomer + i*numAnglesPerMolecule
			if(d[0] == -1):
				tmp.append(d)
		if(j != NUM_MONOMERS-1):
			for t in tmp:
				dataFile.write(str(n+1) + " " + str(t[1]) + " " + str(t[2] + start) + " " + str(t[3] + start) + " " + str(t[4] + start) + "\n")
				n += 1
dataFile.write("\n")

'''
dataFile.write("Ellipsoids\n\n")

for mol in range(NUM_MOLECULES):
	for i in range(NUM_MONOMERS):
		for d in monomer_data_atoms:
			if(d[6] == 1):
				dataFile.write(str(d[0] + mol*numAtomPerChain + i*numAtomsPerMonomer) + " " + str(d[8]) + " " + str(d[9]) + " " + str(d[10]) + " 1 0 0 0\n")
dataFile.write("\n")

dataFile.write("CoreIDs\n\n")
for i in range(NUM_MOLECULES):
	for j in range(NUM_MONOMERS):
		for d in monomer_data_rigid:
			dataFile.write(str(d[0] + i*numAtomPerChain + j*numAtomsPerMonomer) + " " + str(d[1] + i*monomer_disks*NUM_MONOMERS + j*monomer_disks) + "\n")
'''
