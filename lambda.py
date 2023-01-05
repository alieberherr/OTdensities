import numpy as np
import utils
import argparse
import csv
import os

parser = argparse.ArgumentParser(description="Computation of Lambda descriptor")
parser.add_argument('--dir', default=os.getcwd() + '/', help="working directory for orbitals")
parser.add_argument('--csvfile', help='name of CSV file with excitation data')
parser.add_argument('--orbitals', help='home directory of orbitals')
parser.add_argument('--results', help='results directory')
parser.add_argument('--format', help='which program is used to obtain orbitals')

def overlap(phia, phib):
	''' calculates the overlap E(|phi_i| |phi_a|)
	parameters:
	phia: name of the file where the absolute values of the source orbital are stored
	phib: name of the file where the absolute values of the target orbital are stored
	returns:
	the value of the overlap between the moduli of the two densities.
	'''
	grid, vals1, stats1 = utils.cube_to_array(phia)
	_, vals2, stats2 = utils.cube_to_array(phib)
	return utils.integrate3D(grid, np.abs(vals1)*np.abs(vals2))

if __name__ == '__main__':
	args = parser.parse_args()

	# initialise where the results are to be stored
	resultsfile = args.results

	with open(resultsfile, 'w') as wfile:
		writer = csv.DictWriter(wfile, fieldnames=["Molecule", "Excitation", "Functional", "Lambda1", "Lambda2"])
		writer.writeheader()
		with open(args.csvfile,'r') as file:
			reader = csv.DictReader(file)
			# every line corresponds to one excitation
			for line in reader:
				failed=False
				if line["no contr"] == '':
					continue
				molecule = line["\ufeffMolecule"]
				# molecule = line["Molecule"]
				excitation = line["Excitation"]
				functional = line["Functional"]
				if functional == "CAM-B3LYP":
					functional = "CAMB3LYP"
				no_contr = line["no contr"]

				versions = 1
				if "Delta" in excitation or "Sigma-" in excitation:
					versions = 2
				Lambda = np.zeros(versions)

				for i in range(versions):
					for j in range(int(no_contr)):
						phia = line["occ" + str(j + 1)]
						phib = line["virt" + str(j + 1)]
						c = line["contr" + str(j + 1)]
						if args.format=='turbomole':
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


					if args.format=='orca':
						filenamea = args.dir + molecule + "/" + args.orbitals + molecule + "_" + functional + ".mo" + phia + ".cube"
						filenameb = args.dir + molecule + "/" + args.orbitals + molecule + "_" + functional + ".mo" + phib + ".cube"
					elif args.format=='turbomole':
						phia = phia.replace(" ","")
						phib = phib.replace(" ","")
						filenamea = args.dir + molecule + "/" + functional + args.orbitals + phia + ".cub"
						filenameb = args.dir + molecule + "/" + functional + args.orbitals + phib + ".cub"

					if os.path.isfile(filenamea) and os.path.isfile(filenameb):
						ov = overlap(filenamea, filenameb)
						Lambda[i] += ov*float(c)
					else:
						print("Failed:",molecule, excitation, functional, phia, phib, versions)
						print("Looking for files:")
						print(filenamea)
						print(filenameb)
						failed = True


				if not failed:
					for i in range(versions):
						writer.writerow({"Molecule": molecule, "Excitation": excitation, "Functional": functional, f"Lambda%i"%(i+1): Lambda[i]})
