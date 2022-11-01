import numpy as np
import matplotlib.pyplot as plt

'''
ITEM: TIMESTEP
100000
ITEM: NUMBER OF ATOMS
2700
ITEM: BOX BOUNDS pp pp pp
-5.0000000000000000e+02 5.0000000000000000e+02
-5.0000000000000000e+02 5.0000000000000000e+02
-5.0000000000000000e+02 5.0000000000000000e+02
ITEM: ATOMS id type xu yu zu c_op[1] c_op[2] c_op[3] c_op[4] c_shape[1] c_shape[2] c_shape[3]
'''

bondSite_type = 2

ell = open("traj_ell.5.dump", "r")
all_atoms = open("traj_all.5.dump", "r")
combined = open("viz.dump", "w")

# Finding number of ellipsoid particles
ell_line = ""
for i in range(4):
	ell_line = ell.readline()
num_ell = int(ell_line)
ell.seek(0)

# Finding the number of total particles
all_line = ""
for i in range(4):
        all_line = all_atoms.readline()
num_all = int(all_line)
all_atoms.seek(0)

ell_line = ell.readline()
all_line = all_atoms.readline()
cond = True
while cond:
	if ell_line == "" or all_line == "":
		break
	for i in range(1):
		ell_line = ell.readline()
		all_line = all_atoms.readline()
	if int(ell_line) != int(all_line):
		print("error: timesteps do not match up")
	combined.write("ITEM: TIMESTEP\n")
	combined.write(ell_line)
	combined.write("ITEM: NUMBER OF ATOMS\n")
	combined.write(str(num_all) + "\n")
	combined.write("ITEM: BOX BOUNDS pp pp pp\n")
	for i in range(4):
		ell_line = ell.readline()
		all_line = all_atoms.readline()
	for i in range(3):
		combined.write(all_line)
		ell_line = ell.readline()
		all_line = all_atoms.readline()
	combined.write(ell_line)
	#cond = False
	ellipsoids = {}
	all_things = {}
	for i in range(num_ell):
		ell_line = ell.readline()
		ell_split = ell_line.split()
		key = int(ell_split[0])
		data = [float(ell_split[5]), float(ell_split[6]), float(ell_split[7]), float(ell_split[8]), float(ell_split[9]), float(ell_split[10]), float(ell_split[11])]
		ellipsoids[key] = data

	for i in range(num_all):
		all_line = all_atoms.readline()
		all_split = all_line.split()
		key = int(all_split[0])
		data = [float(all_split[1]), float(all_split[2]), float(all_split[3]), float(all_split[4])]
		all_things[key] = data

	for i in range(1, num_all+1):
		if(all_things[i][0] == bondSite_type):
			combined.write(str(i) + " " + 	str(int(all_things[i][0])) + " " + str(all_things[i][1]) + " " +  str(all_things[i][2]) + " " + str(all_things[i][3]) + " 0 0 0 1 0.01 0.01 0.01\n")
		else:
			combined.write(str(i) + " " +  str(int(all_things[i][0])) + " " + str(all_things[i][1]) + " " +  str(all_things[i][2]) + " " + str(all_things[i][3]) + " " + str(ellipsoids[i][0]) + " " + str(ellipsoids[i][1]) + " " + str(ellipsoids[i][2]) + " " + str(ellipsoids[i][3]) + " " + str(ellipsoids[i][4]) + " " + str(ellipsoids[i][5]) + " " + str(ellipsoids[i][6]) + "\n")


	ell_line = ell.readline()
	all_line = all_atoms.readline()
