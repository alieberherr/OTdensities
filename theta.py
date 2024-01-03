import torch
from geomloss import SamplesLoss
import numpy as np
import numpy.linalg
import utils
import argparse
import csv
import os


parser = argparse.ArgumentParser(description="Computation of Theta descriptor")
parser.add_argument('--dir', default=os.getcwd() + '/', help="working directory for orbitals")
parser.add_argument('--csvfile', help='name of CSV file with excitation data')
parser.add_argument('--orbitals', help='home directory of orbitals')
parser.add_argument('--results', help='results directory')
parser.add_argument('--format', help='which program is used to obtain orbitals')


def sinkhorndiv(phia,phib):
	''' Calculates the Sinkhorn divergence of |phi_a|**2 to |phi_b|**2 based on the
	geomloss package.
	'''
	grid,vals1,stats1 = utils.cube_to_array(phia)
	_,vals2,stats2 = utils.cube_to_array(phib)
	vals1 = vals1/np.sqrt(np.sum(vals1**2))
	vals2 = vals2/np.sqrt(np.sum(vals2**2))
	vals1 = vals1**2+1e-16
	vals2 = vals2**2+1e-16
	mxdist = np.linalg.norm(grid[0,:]-grid[-1,:])**2
	eps = mxdist/1e4
	Loss = SamplesLoss("sinkhorn", p=2, blur=eps, scaling=0.8)
	
	v1 = torch.from_numpy(vals1)
	v2 = torch.from_numpy(vals2)
	x = torch.from_numpy(grid)
	y = torch.clone(x)
	return Loss(v1,x,v2,y).item(),eps


if __name__=="__main__":
	use_cuda = torch.cuda.is_available()
	dtype = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor
	
	args = parser.parse_args()

	# initialise where the results are to be stored
	resultsfile = args.results

	with open(resultsfile, 'w') as wfile:
		writer = csv.DictWriter(wfile, fieldnames=["Molecule", "Excitation", "Functional", "Theta1", "Theta2", "eps"])
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
				print("Moleule:",molecule)
				print("Excitation:",excitation)
				print("Functional:",functional)
				print("Contributions:")

				versions = 1
				if "Delta" in excitation or "Sigma-" in excitation or "Sigmau-" in excitation:
					versions = 2
				Theta = np.zeros(versions)

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
							th,eps = sinkhorndiv(filenamea, filenameb)
							print(f"Theta %6.4f with contribution %5.3f"%(th,float(c)))
							Theta[i] += th*float(c)
						else:
							print("Failed:",molecule, excitation, functional, phia, phib, versions)
							print("Looked for files:")
							print(filenamea)
							print(filenameb)
							failed = True
							eps = 1000
							break


				if not failed:
					if (versions==1):
						writer.writerow({"Molecule": molecule, "Excitation": excitation, "Functional": functional, "Theta1": Theta[0], "eps":eps})
					else:
						writer.writerow({"Molecule": molecule, "Excitation": excitation, "Functional": functional, "Theta1": Theta[0], "Theta2": Theta[1], "eps":eps})
