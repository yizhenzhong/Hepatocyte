file = open("replicate.v7.log")
output = open("replicate.process.v7.log","w")
temp = []
for line in file:

	if "GTEx" in line:
		if len(temp) != 0:
			output.write("\t".join(temp) + "\n")
			temp = []
			temp.append(line.strip().split("/")[-1])
		else:
			temp.append(line.strip().split("/")[-1])
		
	elif "temp.txt" in line:
		temp.append(line.split()[0])
	else:	
		
		temp.append(line.split()[0])
		
output.write(" ".join(temp) + "\n")
output.close()
