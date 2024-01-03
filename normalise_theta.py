import numpy as np
import pandas as pd
import argparse
import csv

parser = argparse.ArgumentParser(description="Normalising Theta descriptor by <r^2> of the orbitals (assumes <r^2> has been precomputed)")
parser.add_argument('--thetafile', help='name of CSV file with Theta data')
parser.add_argument('--excfile', help='name of CSV file with excitation data')
parser.add_argument('--orbfile', help='home directory of orbitals')
parser.add_argument('--resfile', help='results file')

if __name__=="__main__":
	args = parser.parse_args()

	orbdat = pd.read_csv(args.orbfile)
	thetadat = pd.read_csv(args.thetafile)

	f = open(args.resfile,'w')
	f.write("Molecules,Functional,Excitation,Type,Theta 1 (orig),Theta 2 (orig),Theta 1 (norm),Theta 2(norm)\n")

	reader = csv.DictReader(open(args.excfile))
	for line in reader:
		molecule = line["\ufeffMolecule"]
		excitation = line["Excitation"]
		functional = line["Functional"]
		extype = line["Type"]

		no_contr = int(line["no contr"])

		versions = 1
		if "Delta" in excitation or "Sigma-" in excitation or "Sigmau-" in excitation:
			versions = 2
		Theta = np.zeros(versions)
		Theta_norm = np.zeros(versions)

		for i in range(versions):
			for j in range(no_contr):
				phia = line["occ%i"%(j+1)]
				phib = line["virt%i"%(j+1)]
				c = float(line["contr%i"%(j+1)])

				if "e" in phia and "e" in phib:
					if i==0:
						phia += "1"
						phib += "1"
					if i==1:
						phia += "1"
						phib += "2"
					if i==2:
						phia += "2"
						phib += "1"
					if i==3:
						phia += "2"
						phib += "2"
				elif "e" in phia:
					if i==0:
						phia += "1"
					if i==1:
						phia += "2"
				elif "e" in phib:
					if i==0:
						phib += "1"
					if i==1:
						phib += "2"

				phia = phia.replace(" ","")
				phib = phib.replace(" ","")

				thetacurr = thetadat[(thetadat["Molecule"]==molecule) & (thetadat["Functional"]==functional) & (thetadat["occ"]==phia) & (thetadat["virt"]==phib)]
				ava = orbdat[(orbdat["Molecule"]==molecule) & (orbdat["Functional"]==functional) & (orbdat["Orbital"]==phia)]
				avb = orbdat[(orbdat["Molecule"]==molecule) & (orbdat["Functional"]==functional) & (orbdat["Orbital"]==phib)]

				if (len(thetacurr)==0):
					continue

				thetacurr = thetacurr["S"].item()
				ava = ava["<r^2>"].item()
				avb = avb["<r^2>"].item()
				Theta[i] += c*thetacurr
				Theta_norm[i] += c*thetacurr/np.sqrt(ava*avb)

		if versions==1:
			f.write("%s,%s,%s,%s,%f,%f,%f,%f\n"%(molecule,functional,excitation,extype,Theta[0],float("nan"),Theta_norm[0],float("nan")))
		if versions==2:
			f.write("%s,%s,%s,%s,%f,%f,%f,%f\n"%(molecule,functional,excitation,extype,Theta[0],Theta[1],Theta_norm[0],Theta[1]))



