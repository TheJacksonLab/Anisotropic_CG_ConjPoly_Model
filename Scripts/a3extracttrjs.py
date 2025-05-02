
start, stop = 0, 51

##### OVITO FILE #####
f_in_name = 'ovito.trj'
f_in = open(f_in_name, 'r')
f_out = open('ovito{}-{}.trj'.format(start,stop),'w')

trjs = 0
line = f_in.readline()
while line != '':
    if 'ITEM: TIMESTEP' in line:
        trjs += 1
        f_in.readline()
        f_in.readline()
        n = int(f_in.readline())
    line = f_in.readline()
f_in.seek(0)
trj = 0
lines = ''
line = f_in.readline()
while line != '':
    if trj >= start and trj < stop:
        lines += line
        for i in range(n+8): #check 9
            #line = f_in.readline()
            lines += f_in.readline()
        f_out.write(lines)
        lines = ''
    else:
        for i in range(n+8): #check 9
            line = f_in.readline()
    trj += 1
    line = f_in.readline()
f_in.close
f_out.close

##### NBBPOTS FILE #####
f_in_name = 'dComp_pe.out'
f_in = open(f_in_name, 'r')
f_out = open('nbbpots{}-{}.out'.format(start,stop),'w')

trjs = 0
line = f_in.readline()
while line != '':
    if 'ITEM: TIMESTEP' in line:
        trjs += 1
        f_in.readline()
        f_in.readline()
        n = int(f_in.readline())
    line = f_in.readline()
f_in.seek(0)
trj = 0
lines = ''
line = f_in.readline()
while line != '':
    if trj >= start and trj < stop:
        lines += line
        for i in range(n+8): #check 9
            #line = f_in.readline()
            lines += f_in.readline()
        f_out.write(lines)
        lines = ''
    else:
        for i in range(n+8): #check 9
            line = f_in.readline()
    line = f_in.readline()
    trj += 1
f_in.close
f_out.close
