#python gene_puller.py input output
from sys import argv
script,inputfile,output = argv
f = open(output,"w")
with open(inputfile, "r") as ins:
	count = 0
	for line in ins:
		stringer = []
		if count == 0:
			f.write(line)
			count = count + 1
		else:
			array = line.split("\t")
			if array[2] == "gene":
				f.write(line)
f.close()
ins.close()        
  
