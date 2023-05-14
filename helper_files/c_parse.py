import sys 
import os
import numpy as np

file_name = str(sys.argv[1])

data = np.loadtxt(file_name).flatten(order="F")

write_file_split= file_name.split("_")
write_file = str(write_file_split[0])+"_"+str(write_file_split[1])+"_C.out"

with open(write_file,'w') as out_file:
		for i in data:
			out_file.write(str("{:.33f}".format(float(i))+' \n'))