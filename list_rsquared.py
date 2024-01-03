import csv
import numpy as np

exfile = "../data_final/turbomole_SOtransitions.csv"

# compile lists of all pairs phi_i, phi_a to calculate S
orbpairs = {}

with open(exfile,'r') as file:
	reader = csv.DictReader(file)
	for j,line in enumerate(reader):
		molecule = line["\ufeffMolecule"]
		functional = line["Functional"]
		no_contr = int(line["no contr"])

		currpairs = []
		for i in range(no_contr):
			phia = line["occ%i"%(i+1)]
			phib = line["virt%i"%(i+1)]
			phia = phia.replace(" ","")
			phib = phib.replace(" ","")
			# if 'e' in orbitals: use all versions
			versions = 1
			if 'e' in phia and 'e' in phib:
				versions = 4
			elif 'e' in phia or 'e' in phib:
				versions = 2

			for i in range(versions):
				phiaf = phia
				phibf = phib
				if "e" in phia and "e" in phib:
					if i==0:
						phiaf += "1"
						phibf += "1"
					if i==1:
						phiaf += "1"
						phibf += "2"
					if i==2:
						phiaf += "2"
						phibf += "1"
					if i==3:
						phiaf += "2"
						phibf += "2"
				elif "e" in phia:
					if i==0:
						phiaf += "1"
					if i==1:
						phiaf += "2"
				elif "e" in phib:
					if i==0:
						phibf += "1"
					if i==1:
						phibf += "2"
				currpairs.append((phiaf,phibf))

		if (molecule,functional) in orbpairs:
			# merge existing list with currpairs
			prevpairs = orbpairs[(molecule,functional)]
			newpairs = list( set(currpairs) | set(prevpairs) )
			orbpairs[(molecule,functional)] = newpairs
		else:
			# initialise
			orbpairs[(molecule,functional)] = currpairs

# write the pairs needed to disk
f = open("../data_final/rsquared/SOT.csv",'w')
f.write("Molecule,Functional,occ,virt\n")
for (molecule,functional) in orbpairs:
	pairs = orbpairs[(molecule,functional)]
	for i in range(len(pairs)):
		f.write("%s,%s,%s,%s\n"%(molecule,functional,pairs[i][0],pairs[i][1]))
f.close()

# compile list of all orbitals
orbs = {}
for key in orbpairs:
	pairs = orbpairs[key]
	orb = [tmp for tmpl in pairs for tmp in tmpl]
	orbs[key] = [tmp for tmp in set(orb)]

# write the orbitals needed to disk
f = open("../data_final/rsquared/orbs.csv",'w')
f.write("Molecule,Functional,orb\n")
for (molecule,functional) in orbs:
	orb = orbs[(molecule,functional)]
	for i in range(len(orb)):
		f.write("%s,%s,%s\n"%(molecule,functional,orb[i]))
f.close()
