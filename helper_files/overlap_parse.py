import os
import sys

file_name = str(sys.argv[1])

write_file_split= file_name.split("_")
write_file = str(write_file_split[0])+"_"+str(write_file_split[1])+"_overlap.out"

with open(file_name,'r') as atom_file:
	with open(write_file,'w') as out_file:
		for line in atom_file:
			res = line.split(',')
			for i in res:
				out_file.write(str("{:.33f}".format(float(i))+' \n'))