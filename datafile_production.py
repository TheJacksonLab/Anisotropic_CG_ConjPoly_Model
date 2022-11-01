# Script for generating a LAMMPS data file from information 
# specified in a .mon file
import numpy as np

#############################################################
#############################################################
##########        VARIABLES THE USER CAN SET       ##########
#############################################################
#############################################################

# Number of polymers to include
# Change to NUM_POLYMERS
NUM_MOLECULES = 1

# Number of monomers to include in each polymer
# If two monomers are defined in the monomer file
# then this number should be halved
NUM_MONOMERS = 32

# Distance a particle in one monomer to the
# corresponding particle in the next monomer
# In other words, the unit cell size
monomerOffset = 2.4

# Distance to separate different polymers by
# This parameter is necessary because the LAMMPS
# Replicate command does not work with rigid bodies defined
moleculeOffset = 15

# File name of the monomer information file
mon_file_name = "v5_nosc.mon"

# File name of the output data file
data_file_name = "v5_nosc.data"

# Dimensions of the box
# Box goes from 0 to box_i
box_x = 100.0
box_y = 100.0
box_z = 100.0

#############################################################
#############################################################
##########     END OF USER-SPECIFIED VARIABLES     ##########
#############################################################
#############################################################

mon_data = open(mon_file_name, "r")

dataFile = open(data_file_name,"w")

# Defining arrays to store information from
# each of the monomer file categories
monomer_data_atoms = []
monomer_data_angles = []
monomer_data_bonds = []
monomer_data_rigid = []

# Defining variables to contain the total number
# of disks, bonds, and angles in each monomer
monomer_disks = 0
connecting_bonds = 0
connecting_angles = 0

# Reading in data from the monomer file and 
# storing in the variables definied above
line = mon_data.readline()
while(line != "Particles:\n"):
	line = mon_data.readline()

line = mon_data.readline()

# Gathering information from the Particles section
while(line != "Angles:\n"):
	splitline = line.split()
	if(int(splitline[6]) == 1):
		monomer_disks += 1
	temp = [int(splitline[0]), int(splitline[1]), float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), int(splitline[6]), float(splitline[7]), float(splitline[8]), float(splitline[9]), float(splitline[10])]
	monomer_data_atoms.append(temp)
	line = mon_data.readline()

line = mon_data.readline()

# Gathering information from the Angles section
while(line != "Dihedrals/Bonds:\n"):
	splitline = line.split()
	if(int(splitline[0]) == -1):
		connecting_angles += 1
	temp = [int(splitline[0]), int(splitline[1]), int(splitline[2]), int(splitline[3]), int(splitline[4])]
	monomer_data_angles.append(temp)
	line = mon_data.readline()

line = mon_data.readline()

# Gathering information from the Dihedrals/Bonds section
while(line != "Rigid:\n"):
	splitline = line.split()
	if(int(splitline[0]) == -1):
		connecting_bonds += 1
	temp = [int(splitline[0]), int(splitline[1]), int(splitline[2]), int(splitline[3])]
	monomer_data_bonds.append(temp)
	line = mon_data.readline()

line = mon_data.readline()

# Gathering information from the Rigid section
while(line != ""):
	splitline = line.split()
	temp = [int(splitline[0]), int(splitline[1])]
	monomer_data_rigid.append(temp)
	line = mon_data.readline()

# Number of ellipsoids in whole simulation
numEllipsoids = NUM_MOLECULES*NUM_MONOMERS*monomer_disks

# Number of monomers in the whole simulation
numMonomers = NUM_MOLECULES*NUM_MONOMERS

# Number of atoms per monomer
numAtomsPerMonomer = len(monomer_data_atoms)

# Number of atoms per chain
numAtomPerChain = NUM_MONOMERS*len(monomer_data_atoms)

# Number of atoms in the whole simulation
numAtom = numAtomPerChain * NUM_MOLECULES

# Number of bonds per monomer
numBondsPerMonomer = len(monomer_data_bonds)

# Number of bonds per molecule
numBondsPerMolecule = NUM_MONOMERS*len(monomer_data_bonds) - connecting_bonds

# Number of bonds in the whole simulation
numBonds = numBondsPerMolecule*NUM_MOLECULES

# Number of angles per monomer
numAnglesPerMonomer = len(monomer_data_angles)

# Number of angles per molecule
numAnglesPerMolecule = NUM_MONOMERS*len(monomer_data_angles) - connecting_angles

# Number of angles in the whole simulation
numAngles = numAnglesPerMolecule*NUM_MOLECULES

# Number of atom, bond, and angle types in the simulation
atomTypes = int(np.amax(np.array(monomer_data_atoms)[:,1]))
bondTypes = int(np.amax(np.array(monomer_data_bonds)[:,1]))
angleTypes = int(np.amax(np.array(monomer_data_angles)[:,1]))

# Writing header of datafile
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

dataFile.write(str(0) + " " + str(box_x) + " " + "xlo xhi\n")
dataFile.write(str(0) + " " + str(box_y) + " " + "ylo yhi\n")
dataFile.write(str(0) + " " + str(box_z) + " " + "zlo zhi\n")
dataFile.write("\n")

# Writing Atoms section of datafile

dataFile.write("Atoms\n\n")

# Loop through each polymer
for mol in range(NUM_MOLECULES):
	# Starting position of the current polymer
	startX = 0.0
	startY = 0.0 + moleculeOffset*mol # Polymers are offset form each other in y-dir
	startZ = 0.0
	# Loop through all monomers in each polymer
	for i in range(NUM_MONOMERS):
		x = startX + monomerOffset*i
		y = startY
		z = startZ
		# Write data in format that corresponds to a atom style of hybrid full ellipsoid
		#ID type x y z moleculeID charge ellipsoidFLAG density
		for d in monomer_data_atoms:
			dataFile.write(str(d[0] + mol*numAtomPerChain + i*numAtomsPerMonomer) + " " + str(d[1]) + " " + str(d[2]+x) + " " + str(d[3]+y) + " " + str(d[4]+z) + " " + str(mol+1) + " " + str(d[5]) + " " + str(d[6]) + " " + str(d[7]) + "\n")


dataFile.write("\n")

# Writing Bonds section of datafile

dataFile.write("Bonds\n\n")

# Loop through each polymer
for i in range(NUM_MOLECULES):
	# Loop through each monomer
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

# Writing Angles section of datafile

dataFile.write("Angles\n\n")

# Loop through each polymer
for i in range(NUM_MOLECULES):
	# Loop through each monomer
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

# Writing Ellipsoids section of datafile

dataFile.write("Ellipsoids\n\n")


# Loop through each polymer
for mol in range(NUM_MOLECULES):
	# Loop through each monomer
	for i in range(NUM_MONOMERS):
		for d in monomer_data_atoms:
			if(d[6] == 1):
				dataFile.write(str(d[0] + mol*numAtomPerChain + i*numAtomsPerMonomer) + " " + str(d[8]) + " " + str(d[9]) + " " + str(d[10]) + " 1 0 0 0\n")
dataFile.write("\n")

# Writing the CoreIDs section of datafile

dataFile.write("CoreIDs\n\n")

# Loop through each polymer
for i in range(NUM_MOLECULES):
	# Loop through each monomer
	for j in range(NUM_MONOMERS):
		for d in monomer_data_rigid:
			dataFile.write(str(d[0] + i*numAtomPerChain + j*numAtomsPerMonomer) + " " + str(d[1] + i*monomer_disks*NUM_MONOMERS + j*monomer_disks) + "\n")
