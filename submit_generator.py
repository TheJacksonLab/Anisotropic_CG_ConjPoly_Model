#script for submitting multiple jobs

import subprocess
import numpy as np
from numpy import random

#order: [GB Well depth, Angle, Dihedral, random number]
#you can change the values in this list to run different jobs

number_jobs = 8
starts = [10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000, 19000000]
random_seeds = random.randint(9998, size=(number_jobs,3))+1
depths = np.array([1.0])
angles = np.array([1.25])
dihedrals = np.array([3.0])
job_numbers = np.linspace(1,number_jobs,number_jobs).astype(int)
parameters = np.array([depths, angles, dihedrals, job_numbers])
#empty list for storing the names of the directories
names = []

#MAKE DIRECTORIES FOR EACH SET OF PARAMETERS AND FILL AND ARRAY OF NAMES OF EACH FOR EASY ACCESS
for i in range(len(parameters[0])):
	for j in range(len(parameters[1])):
		for k in range(len(parameters[2])):
			for l in range(len(parameters[3])):
				names.append( "W_" + str(float(parameters[0][i])) + "_A_" + str(float(parameters[1][j])) + "_D_" + str(float(parameters[2][k])) + "_J_" + str(int(parameters[3][l])) + "_v5")
				subprocess.call(["mkdir", names[-1]])

#FUNCTION FOR WRITING THE IN.CONINTUE FILES
def write_continue(file, name):
	splitname = name.split("_")
	job = int(splitname[7]) - 1
	file.write("read_restart ../v5_64/restarts_v5/equil." + str(starts[job]) + ".rst\n\n")
	file.write("pair_style gayberne 1.0 1.0 1.0 3.0\n")
	file.write("fix prop all property/atom i_coreid\n")
	file.write("special_bonds lj 0.0 1.0 1.0\n\n")
	file.write("pair_coeff 1 1 " + splitname[1] + " 1.0 0.1 1.0 0.25 0.1 1.0 0.25\n")
	file.write("pair_coeff 2 2 0.0 1.0 1 1 1 1 1 1\n")
	file.write("pair_coeff * 2 0.0 1.0 0 0 0 0 0 0\n\n")
	file.write("bond_coeff 1 harmonic 0.0 0.5\n")
	file.write("bond_coeff 2 harmonic 200.0 0.2\n")
	file.write("bond_coeff 3 eldihedral 0 " + splitname[5] + " 0 0\n\n")
	#file.write("bond_coeff 3 eldihedral 0 0 0 " + splitname[5] + "\n\n")
	file.write("angle_coeff 1 " + splitname[3] + " 180.0\n\n")
	file.write("group disks type 1 2\ngroup ell type 1\ngroup bondSites type 2\n\n")
	#file.write("neigh_modify delay 0 every 1 check no\n")
	##NEw
	file.write("neigh_modify delay 1 exclude group ell bondSites\n")
	#file.write("delete_bonds all bond 1\n\n")
	###
	file.write("variable temps world .085 .1 .115 .13 .145 .16 .175 .19 .205 .22 .235 .25 .265 .28 .3 .33 .36 .39 .42 .45 .48 .51 .54 .57 .6 .65 .7 .75 .8 .85 .9\n")
	file.write("variable init_temps world 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2\n\n")
	file.write("fix 1 all rigid/nve/small custom i_coreid\n")
	# TESTING
	file.write("fix 2 ell langevin ${init_temps} ${temps} 2 " + str(random_seeds[job][0]) + " angmom 1.0\n")
	file.write("fix drift all momentum 10000 linear 1 1 1 angular rescale\n\n")
	file.write("compute ta ell temp/asphere dof all\n")
	file.write("compute rot all erotate/rigid 1\n")
	file.write("compute kin all ke/rigid 1\n")
	file.write("compute rg ell gyration\n")
	file.write("compute myRDF ell rdf 150\n")
	file.write("compute op ell property/atom quati quatj quatk quatw\n")
	file.write("compute shape ell property/atom shapex shapey shapez\n")
	file.write("variable rep world 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30\n")
	file.write("fix outrg ell ave/time 1 1 25000 c_rg file gyr.${rep}.txt\n")
	file.write("fix outrdf ell ave/time 1 1 25000 c_myRDF[*] file rdf.${rep}.txt mode vector\n")
	file.write("dump 1 bondSites custom 25000 traj_bs.${rep}.dump id type xu yu zu\n\n") #OUTPUTS A FILE NEEDED FOR ORDER PARAMETER CALCULATION
	file.write("dump 2 ell custom 25000 traj_ell.${rep}.dump id type xu yu zu c_op[1] c_op[2] c_op[3] c_op[4] c_shape[1] c_shape[2] c_shape[3]\n")
	file.write("dump 3 all custom 25000 traj_all.${rep}.dump id type xu yu zu\n")
	file.write("fix averages all ave/time 10 50 1000 c_thermo_temp c_ta c_thermo_pe c_kin file thermo.${rep}.avg\n\n")
	file.write("timestep 0.00075\n\nthermo_style custom step temp c_ta epair ebond eangle etotal press pe ke c_rot c_kin\n\nthermo 10000\n\n")
	file.write("run 1000000\n")
	file.write("unfix 2\n")
	file.write("timestep 0.00075\n")
	# TESTING
	file.write("fix 2 ell langevin ${temps} ${temps} 2 " + str(random_seeds[job][0]) + " angmom 1.0\n\n")
	file.write("temper 100000000 100000 ${temps} 2 " + str(random_seeds[job][1]) + " " + str(random_seeds[job][2]) + "\n\n")

#FUNCTION FOR WRITING THE SUBMIT SCRIPTS
def write_submit(file, name):
	file.write("#!/bin/bash\n\n")
	file.write("#SBATCH --mail-user=aecohen@uchicago.edu\n") #YOU WILL NEED TO CHANGE THIS LINE
	file.write("#SBATCH --mail-type=ALL\n")
	file.write("#SBATCH --workdir=/scratch/midway2/aecohen/gb_chain/" + name + "\n") #YOU WILL NEED TO CHANGE THIS LINE
	file.write("#SBATCH --nodes=2\n")
	file.write("#SBATCH --ntasks-per-node=16\n")
	file.write("#SBATCH --exclusive\n")
	file.write("#SBATCH --job-name=v5_A25\n")
	file.write("#SBATCH --time=5-00:00:00\n")
	file.write("#SBATCH --partition=depablo-sandyb\n")
	file.write("#SBATCH --qos=normal\n\n")
	file.write("module unload openmpi\nmodule load openmpi/2.0.1\n")
	file.write("mpirun -np 31 ../lmp_R5 -partition 31x1 -in in.continue > savefile.out\n")
	file.write("wait")

#SUBMITTING ALL OF THE SCRIPTS
for i in names:
	temp_cont = open("./" + i + "/in.continue", "w")
	temp_submit = open("./" + i + "/submit.sbatch", "w")
	write_continue(temp_cont, i)
	write_submit(temp_submit, i)
	temp_cont.close()
	temp_submit.close()
	subprocess.call(["sbatch", "./" + i + "/submit.sbatch"]) 	

