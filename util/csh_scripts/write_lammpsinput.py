from argparse import ArgumentParser
from math import sqrt,floor

def parse_args():
    parser = ArgumentParser(
        description = 'Create input lammps file from pot file')
    parser.add_argument(
        '-i',type=str,action="store",default='in.nafepo4',
        help = 'Input lammps file')
    parser.add_argument(
        '-f',type=str,action="store",default='pot_all.in',
        help = 'Input potfile file')
    parser.add_argument(
        '-n',type=int,action="store",default=6,
         help = 'Number of atoms')
    parser.add_argument(
        '-a',type=int,action="store",default=1,
        help = 'are there angular interactions? 1-yes 0-no')
    parser.add_argument(
        '-o',type=str,action="store",default='in.nafepo4_potfit',
        help = 'Write to this file')
    args = parser.parse_args()
    return args

def isBlank (myString):
    if myString and myString.strip():
        return False
    return True

def row_index(i, M):
    ii = M*(M+1)/2-1-i
    K = floor((sqrt(8*ii+1)-1)/2)
    return M-K

def column_index(i, M):
    ii = M*(M+1)/2-1-i
    K = floor((sqrt(8*ii+1)-1)/2)
    jj = ii - K*(K+1)/2
    return M-jj


#----- MAIN -----#

args = parse_args()
lfile = args.i
pfile = args.f
natoms = args.n
angle = args.a
output = args.o

# read potfit potential file

c = 0
num_int = -1
block = 0
param = []
x = []
with open(pfile) as p:
  #count number of interactions
  for line in p:
     if '#F ' in line:
        data = line.split()
        num_int = float(data[2])
     if '#' in line:
        continue    
     if 'type' in line:
        c +=1
        block = 1
        continue
     if block == 1:
        if isBlank(line):       
	   block  = 0
           param.append(x)        
	   x = []
        elif 'cutoff' in line: 
           continue
        else: 
          data = line.split()
          x.append(float(data[1]))

#Eliminate zero interactions

matter = []
ang = []
cs= []
ic = []
ir = []
vdw = 0
csi = 0
angi = 0
for i in range(len(param)):
    if abs(param[i][0])>0.0:	
       if len(param[i])>2:
           ir.append(int(row_index(i,natoms)))
           ic.append(int(column_index(i,natoms)))
           matter.append(param[i])          	
           vdw += 1
       elif len(param[i])==2 and i< num_int-natoms+1:        
           ir.append(int(row_index(i,natoms)))
           ic.append(int(column_index(i,natoms)))
           cs.append(param[i])	   
           csi += 1     	
       else: 
           angi += 1
           ang.append(param[i])    
 

# Write potential values on template

string_vdw =  'pair_style   born/coul/dsf/cs 0.2  8.0  12.0    # A, rho, sigma=0, C, D\n'+ \
        'pair_coeff   *     *       0.0        1.0       0.0  0.0    0.0\n'

for i in range(vdw):
   string_vdw = string_vdw + 'pair_coeff  {0}   {1}   {2}   {3}   0.0  {4}   0.0\n'.format(ir[i],ic[i],matter[i][0],matter[i][1],matter[i][2]*matter[i][1]**6)

str_bond = 'bond_style harmonic\n'

for i in range(csi):
    str_bond = str_bond + 'bond_coeff {0}  {1}  {2}\n'.format(i+1,cs[i][0],cs[i][1]) 

str_angle = 'angle_style harmonic\n'
if angle>0:
	for i in range(angi):
    	    str_angle = str_angle + 'angle_coeff  {0}  {1}  {2}\n'.format(i+1,ang[i][0],ang[i][1]*57.29577951) 

string = ''

ip=0
ib=0
ia=0
with open(lfile) as l:
  for line in l:
      if 'pair_style' in line:
          ip = 1
          string = string + string_vdw + '\n'
      if 'pair_coeff' in line:
          continue
      if ip == 1 and isBlank(line):
          ip = 0
      if 'bond_style' in line:
          ib = 1
          string = string + str_bond  + '\n'
      if 'bond_coeff' in line:
          continue
      if ib == 1 and isBlank(line):
          ib = 0
      if 'angle_style' in line:
          ia = 1
          if angle == 1:
               if len(ang) == 0:
               	  string = string + str_angle + 'angle_coeff  1  0.0  109.5\n'
	       else:
                  string = string + str_angle + '\n' 
      if 'angle_coeff' in line:
          continue
      if ia == 1 and isBlank(line):
          ia = 0
      if (ia == 0 and ib == 0 and ip == 0): 
          string = string + line  

f = open(output,'w')
f.write(string)
