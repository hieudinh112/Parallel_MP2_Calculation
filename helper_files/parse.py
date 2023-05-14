import os
import sys


file_name = str(sys.argv[1])

write_file_split= file_name[:-4]#.split("_")
write_file = write_file_split+"_MOenergy.out"

grab_lines = False
with open(file_name,'r') as atom_file:
	energy_line = []
	for line in atom_file:
		if line.startswith(" Nuclear Repulsion Energy"):
			grab_lines = True
			continue
		elif line.startswith(' Total QAlloc Memory Limit'):
			grab_lines = False
			continue
		elif line.startswith(' There are       ') and grab_lines:
			res = [int(i) for i in line.split(' ') if i.isdigit()]
			with open(write_file,'w') as out_file:
				#out_file.write('alpha and beta\n')
				#for i in res:
				out_file.write(str(res[0])+"\n")
		elif line.startswith(' SCF   energy '):
			res = line.strip().split(" ")
			with open(write_file,'a') as out_file:
				#out_file.write('\nscf energy\n')
				out_file.write(str(res[-1])+'\n')
		elif line.startswith(" Requested basis set"):
			grab_lines = True
			continue
		elif line.startswith(" A cutoff of  "):
			grab_lines = False
			continue
		elif line.startswith(' There are') and grab_lines:
			res = [int(i) for i in line.split(' ') if i.isdigit()]
			with open(write_file,'a') as out_file:
				#out_file.write('\nnbasis\n')
				out_file.write(str(res[-1])+"\n")


energies = []
grab_lines = False
with open(file_name,'r') as atom_file:
    energy_line = []
    for line in atom_file:
        if line.startswith(' -- Occupied --'):
            grab_lines = True
            continue
        elif line.startswith(' -- Virtual --'):
        	continue
        elif line.startswith(' --------------------------------------------------------------'):
            grab_lines = False
            # if energy_line:
            #     #just checks that we aren't appending an empty list.
            #     energies.append(energy_line)
        if grab_lines: #in python 'is True' is implicit for many types.
            #energy_line.append(line)
            with open(write_file,'a') as out_file:
            	res = line.strip().split(" ")
            	for i in res:
            		if i:
            			out_file.write(str(i)+"\n")


# with open('energies'+file_name,'a') as out_file:
# 	#out_file.write('orbital energies \n')
# 	for energy in energies:
# 		out_file.write(''.join(energy))
# 		break


