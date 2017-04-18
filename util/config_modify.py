import sys
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(
        description = 'Modify existing configs.force files. Operations are performed in the following order:\n \
                       1 - Rename\n   \
                       2 - Duplicate\n   \
                       3 - Add element\n \
		       5 - Reorder\n   \
                       6 - Make forces of ith element zero\n   \
                       7 - Delete species\n   \
                       8 - Split configurations\n   \
                       Using these commands multiple times many different operations are possible')                          
    parser.add_argument(
        '-f',type=str,action="store",default = 'configs.force',
        help = 'name of input configs.force file')
    parser.add_argument(
        '-o',type=str,action="store",default = '',
        help = ' "X Z W Y" Reorder elements W X Y Z in the specified order')
    parser.add_argument(
        '-a',type=str,action="store",default = '',
        help = 'Adds an element in the #C line (does not modify anything else)')
    parser.add_argument(
        '-r',type=str,action="store",default = '',
        help = ' "X Y" Renames element X as element Y')
    parser.add_argument(
        '-x',type=str,action="store",default = '',
        help = '"X Y" Duplicates element X and labels the copy as Y')
    parser.add_argument(
        '-s',type=bool,action="store",default=0,
        help = '1: Splits config file into individual CONTCAR files (maybe useful to perform certain modifications)')
    parser.add_argument(
        '-z',type=str,action="store",default='',
        help = ' "X" Make forces on atom X zero')
    parser.add_argument(
        '-d',type=str,action="store",default='',
        help = ' "X" Delete species X')
    args = parser.parse_args()
    return args



def readConf(fname):
   atoms = []
   with open(fname) as f:
      flg = 0
      config = 0
      for line in f:
         if '#N ' in line:
	     if config > 0:
                conf.append(els)
	        if data0 < nlabels - 1:
                   data0 = nlabels - 1
                   conf.append(data0)
                atoms.append(conf)
             els = 0 
             conf = []
             iflg = 0
             config += 1
	     data = line.split()
             conf.append(config)
             conf.append(int(data[1]))
             continue
         if '#C ' in line and flg == 0:
                labels = line.split()
                labels = labels[1:]
                nlabels = len(labels)
                flg = 1
         if '#F' in line:
              iflg = 1
              continue
         if isBlank (line) and iflg == 1:
             continue
	 if iflg == 1:
             els += 1
             data = line.split()
	     if els == 1:
		data0 = int(data[0])
                conf.append(data0)
		conf.append(1)
	        continue  
             if int(data[0]) != data0:
                conf.append(els-1)
                conf.append(int(data[0]))
                conf.append(els)
                data0 = int(data[0])
	
      conf.append(els)
      if data0 < nlabels - 1:
         data0 = nlabels - 1
         conf.append(data0)      
      atoms.append(conf)
      
      return config,atoms,labels,nlabels
     


def isBlank (myString):
    if myString and myString.strip():
        return False
    return True


###### MAIN #######
 
args = parse_args()
fname = args.f
order = args.o
rname = args.r
dupl = args.x
split = args.s
zeroed = args.z
add = args.a
delete = args.d

delta = 0.0

## Rename
if isBlank(rname)==0:
   string = ''
   print('Renaming')
   species = rname.split()
   if len(species)==1:
      sys.exit('Error renaming: there should be two strings: current atom label and new atom label')
   
   with open(fname) as f:
        for line in f:
            if '#C ' in line:
               if species[0] in line:
                  data = line.split()
                  for i in range(len(data)):
                      if species[0] in data[i]:
                         data[i] = species[1] 
                      string = string + data[i] + ' '
                  string = string + '\n'
                  continue
               else:
                  sys.exit('Error renaming: The name of the element to be renamed is incorrect')
            string = string + line
   g = open(fname,'w')
   g.write(string)
   g.close()

## Duplicate
if isBlank(dupl)==0:
   string = ''
   print('Duplicating')
   str2 = ''
   config,atoms,labels,nlabels = readConf(fname)
   i0 = -1
   for i in range(len(labels)):
       if labels[i] in dupl:
	  i0 = i + 1
          break
   if i0 == -1:
          sys.exit('Error duplicating: The name of the element to be duplicated is incorrect')
   species = dupl.split()
   if len(species)==1:
      sys.exit('Error duplicating: there should be two strings: atom label to be duplicated and label of copy')
   conf = 0
   with open(fname) as f:
        for line in f: 
            if '#N ' in line:     
                iflg = 0
                if conf > 0:
                   string = string + newlines
	        conf += 1
                new = str(len(labels))
                newlines = ''
                data = line.split()
                for j in range(config):
                    if atoms[j][0] == conf:
                       thisconf = atoms[j]
		       i1 = 1 + 3*i0
                       extra = thisconf[i1] - thisconf[i1 - 1] + 1
 		       string = string + '#N ' + str(extra + thisconf[1]) + ' ' + str(data[1]) + '\n'
         	       continue
                continue   
            if '#C ' in line:
               if species[0] in line:
		  string = string + line[0:-1] + ' ' + species[1] + '\n'
                  continue
               else:
                  sys.exit('Error duplicating: The name of the element to be duplicated is incorrect')            
            if '#F' in line:
               iflg = 1
               string = string + line
	       continue
            if isBlank(line) and iflg == 1:
	       continue		   
            string = string + line
            if iflg == 1:
               data = line.split()
               if i0-1 == int(data[0]):
                   newlines = newlines + new + ' ' + line[1:-1] + '\n'       
        string = string + newlines
   g = open(fname,'w')
   g.write(string)
   g.close()


## Add
if isBlank(add)==0:
   string = ''
   print('Adding')
   string = ''
   with open(fname) as f:
        for line in f:
            if '#C ' in line:
               string = string + line[0:-1] + ' ' + add + '\n'               
               continue
            string = string + line
   g = open(fname,'w')
   g.write(string)
   g.close()


## Reorder
if isBlank(order)==0:
   string = ''
   new_ord = []
   print('Reordering')
   config,atoms,labels,nlabels = readConf(fname)
   data = order.split()
   for s in data:
       c = -1
       for q in labels: 
           c += 1
           if s == q:
              new_ord.append(c)		
   with open(fname) as f:
        conf = 0
        iflg = 0
        current_line = 0
        copy_all = [] 
        nextline = []
        idx = []
        for line in f:
            if '#C ' in line:
                data = line.split()
                data = data[1:]         
                data = [data[i] for i in new_ord]
                linec = ''
                for s in data:
                    linec = linec + ' ' + s
                linec = '#C' + linec + '\n'
                copy_all.append(linec)
                current_line = current_line + 1
                nextline.append(current_line)
                continue
            copy_all.append(line)
            if '#N' in line:      
               iflg = 0
               conf += 1
            if '#F' in line:
               current_line = current_line + 1
	       nextline.append(current_line)
               iflg = 1
               continue
            if isBlank(line) and iflg == 1:
               current_line = current_line + 1
               nextline.append(current_line)
               continue
            if iflg == 1:
               this_conf = atoms[conf - 1]
               natoms = this_conf[1]
               for i in range(c+1):
		   for j in range(c+1):
		      i0 = 3*(j+1) - 1
                      if len(this_conf) == i0 + 1:
                          continue
                      if this_conf[i0] == new_ord[i]:
                          nat = this_conf[i0+2] - this_conf[i0+1] + 1
                          nextline.extend(range(current_line+this_conf[i0+1], current_line + this_conf[i0+1] + nat))
                          [idx.append(i) for k in range(0,nat)]                          
               current_line = current_line + natoms
               iflg  = 2
               continue
            if iflg == 0:
               current_line = current_line + 1                    
               nextline.append(current_line)
               continue
   current_line = -1
   iflg = 0
   for i in range(len(nextline)): 
        i0 = nextline[i] - 1
        line  = copy_all[i0]
        if '#N ' in line:
           iflg = 0
           string = string + line
           continue
        if '#F' in line:
           iflg = 1
           string = string + line
           continue
        if isBlank(line) and iflg == 1:
            string = string + line
            continue
        if iflg == 1:
           current_line += 1
	   string = string + str(idx[current_line]) + ' '+ line[1:] 
           continue
        string = string + line

   f = open(fname,'w')
   f.write(string)
   f.close()

## Zero
if isBlank(zeroed)==0:
   string = ''
   print('Zeroing')
   config,atoms,labels,nlabels = readConf(fname)
   for i in range(len(labels)):
       if labels[i] == zeroed:
	  i0 = i
   with open(fname) as f:
        for line in f:
  	    if line[0] == str(i0):
	       data = line.split()
               string  = string + data[0] +'     '+data[1]+'     '+data[2]+'     '+data[3]+ \
                         '     0.000000    0.000000    0.000000\n'
               continue
            string = string + line
   f = open(fname,'w')
   f.write(string)
   f.close()

## Deleting
if isBlank(delete)==0:
   string = ''
   print('Deleting')
   config,atoms,labels,nlabels = readConf(fname)
   for i in range(len(labels)):
       if labels[i] == delete:
          i0 = i
          break
   conf = 0
   with open(fname) as f:
        for line in f:
            if '#N' in line:
                conf += 1
                data = line.split()
                this_conf = atoms[conf - 1]
                nat = this_conf[3*i0+4] - this_conf[3*i0+3] + 1
	        string = string + '#N' + '  ' + str(int(data[1]) - nat) + '  ' + data[2] + '\n'
                continue
            if '#C' in line:
		data = line.split()
                data = data[:i0+1] + data[i0+2 :]
                string  = string +  ' '.join(data) + '\n'
                continue
            if line[0] == str(i0):
               continue
            if line[0] > str(i0):
               data = line.split()
               string  = string + str(int(data[0])-1) +'     '+data[1]+ \
                '     '+data[2]+'     '+data[3] +'     '+data[4]+'     '+data[5]+'     '+data[6] + '\n'
               continue
            string = string + line
   f = open(fname,'w')
   f.write(string)
   f.close()

#Split
if split==1:
   string = ''
   print('splitting')
   with open(fname) as f:
        configs = 0
        for line in f:
	    if '#N' in line:
                if configs > 0:
                   g.close()
                new_name = 'CONTCAR_SPLIT_'+str(configs) 
                g = open(new_name,'w')
                configs += 1
            g.write(line)    
  














