import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Normalising Theta descriptor by <r^2> of the orbitals (assumes <r^2> has been precomputed)")
parser.add_argument('--thetafile', help='name of CSV file with Theta data')
parser.add_argument('--excfile', help='name of CSV file with excitation data')
parser.add_argument('--orbfile', help='home directory of orbitals')
parser.add_argument('--resfile', help='results file')

if __name__=="__main__":
	args = parser.parse_args()
	
	orbdat = pd.read_csv(args.orbfile)
	thetadat = pd.read_csv(args.thetafile)
	excdat = pd.read_csv(args.excfile)

	f = open(args.resfile,'w')
	f.write("Molecules,Functional,Excitation,Type,Theta (orig),Theta (norm)\n")

	for (idx,row) in excdat.iterrows():
		molecule = row.loc["Molecule"]
		excitation = row.loc["Excitation"]
		functional = row.loc["Functional"]
		extype = row.loc["Type"]

		ncontr = row.loc["no contr"]
		if ncontr > 1:
			continue
		# find Theta
		theta = thetadat[thetadat["Molecule"] == molecule]
		theta = theta[theta["Functional"] == functional]

		# find overlaps of involved orbitals
		orbi = row.loc["occ1"]
		orbf = row.loc["virt1"]
		orbi = orbi.replace(" ","")
		orbf = orbf.replace(" ","")
		nversions = 1
		if "e" in orbi or "e" in orbf:
			nversions = 2

		for i in range(nversions):
			phia = orbi
			phib = orbf
			if "e" in phia and "e" in phib:
				if i==0:
					phia += "1"
					phib += "1"
				if i==1:
					phia += "1"
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
		
			phia = phia + ".cub"
			phib = phib + ".cub"

			# find overlaps of involved orbitals
			overlap = orbdat[orbdat["Molecule"] == molecule]
			overlap = overlap[overlap["Functional"] == functional]
			ovi = overlap[overlap["Orbital"] == phia]["<r^2>"]
			ovf = overlap[overlap["Orbital"] == phib]["<r^2>"]
			if len(ovi) == 0 or len(ovf) == 0:
				if len(ovi) == 0:
					print("orbi:",orbi)
					print("not found:",phia)
				if len(ovf) == 0:
					print("orbf:",orbf)
					print("not found:",phib)
				continue
			ovi = ovi.item()
			ovf = ovf.item()
	
			# choose the right theta
			if i==0:
				thetacurr = theta[theta["Excitation"] == excitation]["Theta 1 (calc, 0.8)"]
			if i==1:
				thetacurr = theta[theta["Excitation"] == excitation]["Theta 2 (calc, 0.8)"]
			thetacurr = thetacurr.item()
			if np.isnan(thetacurr):
				continue
			# normalise theta
			normtheta = thetacurr/np.sqrt(ovi*ovf)

			f.write(molecule + "," + functional + "," + excitation + "," + extype + ",%f,%f\n"%(thetacurr,normtheta))

	f.close()
