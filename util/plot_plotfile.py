#!/usr/bin/env python3

def parse_args():
    parser = ArgumentParser(
        description = 'Plot functions dumped in potfit potential file')
    parser.add_argument(
        'plotfile',type=str,action="store",
        help = 'plot file written by potfit')
    parser.add_argument(
        'potentialfile',type=str,action="store",
        help = 'potential file used in the same run (input or output)')
    parser.add_argument(
        '-f',type=str,action="store",
        choices=['pdf','svg','png'],
        help = 'format for saving plot figures')
    args = parser.parse_args()
    return args

###### some useful functions
def is_number(s):
    try:
       float(s)
       return True
    except ValueError:
       return False

def notcut(s):
    if 'cutoff' in s:
       return False
    else:
       return True

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

########## parse potential file header
def parseheader(f):
  species = []
  intspe = []
  with open(pfile) as p:
    line = p.readline()
    while line[0] == '#' and line[1] != 'E':
      if line[1] == 'F':
        pformat = int(line.split()[1])
        if pformat == 0:
          num_int = int(line.split()[2])
          print("Analytic potentials format detected with %d interactions" % num_int)
      if line[1] == 'C':
        buf = line.split() 
        for i in range(len(buf)-1):
          species.append(buf[i+1])
      ##try to get interactions labels as written in potential output
      if line[1] == '#':
        buf = line.split() 
        for i in range(len(buf)-1):
          intspe.append(buf[i+1])
        if (len(intspe) == num_int):
          print("Labels for interactions found")
      line = p.readline()

  return num_int,pformat,species,intspe

######## get interactions types and their parameters for analytic potentials
def getinteractions(f):
  block = 0
  param = []
  x = []
  intkind = []
  with open(pfile) as p:
    #count number of interactions
    for line in p:
      if 'type' in line:
        func = line.split()[1]
        intkind.append(func)
        block = 1
        first = True
        continue
      while block == 1:
        if (first):
           x = []
           line2 = line
        else:
           line2 = p.readline()
        first = False
        data = line2.split()
        if len(data)<1:
            block = 0
        elif is_number(data[1]) and notcut(line2):
            x.append(float(data[1]))
      if bool(x):
        param.append(x)

  return intkind,param


#----- MAIN -----#

import sys
import numpy as np
from os.path import isfile
from argparse import ArgumentParser
from math import sqrt,floor

args = parse_args()
ifile = args.plotfile
pfile = args.potentialfile
oformat = args.f

###### check for good args
if ( not isfile(ifile)):
  sys.stderr.write('Plot file %s not found \n' % ifile)
  sys.exit()

if ( not isfile(pfile)):
  sys.stderr.write('Potential file %s not found\n' % pfile)
  sys.exit()

###### read potfit potential file
num_int,pformat,species,intspe = parseheader(pfile)
if pformat == 0:
    intkind, param = getinteractions(pfile)

###### check consistency in previous reads
if (len(intkind) != num_int):
    sys.stderr.write("ERROR: \t%d functions read in header but %d ones 'type' found in file \n" % (num_int,len(intkind)))
    sys.stderr.write("\tUse the -s and/or -f command line switches.\n")
    sys.exit()
elif (len(intkind) != len(param) ):
    sys.stderr.write("ERROR: \tDifferent numbers of functions types and parameters set found\n")
    sys.exit()

for i in range(num_int):
    print(intspe[i],intkind[i],param[i])

sys.exit()

###### read plotfile


#Eliminate zero interactions

matter = []
ang = []
cs = []
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
       elif intkind[i] == 'angharmonic': 
          angi += 1
          ang.append(param[i])    
 
# write potential values on template

string_vdw =  'pair_style   born/coul/dsf/cs 0.2  8.0  12.0    # A, rho, sigma=0, C, D\n'+ \
        'pair_coeff   *     *       0.0        1.0       0.0  0.0    0.0\n'

for i in range(vdw):
  string_vdw = string_vdw + 'pair_coeff  {0}   {1}   {2}   {3}   0.0  {4}   0.0\n'.format(ir[i],ic[i],matter[i][0],matter[i][1],matter[i][2]*matter[i][1]**6)

str_bond = 'bond_style harmonic\n'

for i in range(csi):
  str_bond = str_bond + 'bond_coeff {0}  {1}  {2}\n'.format(i+1,cs[i][0],cs[i][1]) 

if angi > 0:
  str_angle = 'angle_style harmonic\n'
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

