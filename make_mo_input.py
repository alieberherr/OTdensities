import numpy as np
import utils
import os
import csv
import argparse

parser = argparse.ArgumentParser('generate MO plots')
parser.add_argument('--dir', type=str, help='directory (default set to OTproject scratch directory)', default='/cluster/scratch/anninal/OTproject/orbitals/')
parser.add_argument('--ngridpoints', type=str, nargs='+', help='number of grid points', default=[40,40,40])
parser.add_argument('--minpoints', type=str, nargs='+', help='lower limits of grid', default=[])
parser.add_argument('--maxpoints', type=str, nargs='+', help='upper limits of grid', default=[])
parser.add_argument('--molecule', type=str, help='molecule', default='')
parser.add_argument('--functional', type=str, help='functional', default='')
parser.add_argument('--inter', type=int, help='produce input for interactive mode', default=0)

def make_input(mos, molecule, functional, npoints, minpoints, maxpoints):
	out = open(args.dir + molecule + "/" + molecule + "_" + functional + ".inp", "w")
	out.write("! NoIter\n\n%plots\n Format Gaussian_Cube\n")
	if len(minpoints) > 0:
		out.write(" min1 %s\n min2 %s\n min3 %s\n"%(minpoints[0], minpoints[1], minpoints[2]))
	if len(maxpoints) > 0:
		out.write(" max1 %s\n max2 %s\n max3 %s\n"%(maxpoints[0], maxpoints[1], maxpoints[2]))
	for i in range(3):
		out.write(" dim%i %s\n"%(i+1,npoints[i]))
	for mo in mos:
		out.write(" MO(\"formaldehyde_%s.mo%sa.cube\",%i,0);\n"%(functional,mo,int(mo)))
	out.write("end\n\n* xyzfile 0 1 " + molecule + ".xyz")
	out.close()

def make_input_interactive(mos, molecule, functional, npoints):
	out = open(args.dir + molecule + "/plot_" + molecule + "_" + functional + ".txt", "w")
	out.write("5\\n7\\n")
	for i in range(mos.size):
		out.write("2\\n"+mos[i]+"\\n")
		out.write("4\\n"+npoints[0]+" "+npoints[1]+" "+npoints[2]+"\\n")
		out.write("10\\n")
	out.write("11\\n")
	out.close()

if __name__ == '__main__':
	args = parser.parse_args()
	molecule = ''
	excitation = ''
	functional = ''

#	with open('/cluster/scratch/anninal/OTproject/orbitals/ex_refBasis.csv','r') as file:
	with open('/cluster/scratch/anninal/OTproject/orbitals/ex_remaints.csv','r') as file:
		reader = csv.DictReader(file)
		# every line corresponds to one excitation
		mos = []
		for line in reader:
			molecule = line["\xef\xbb\xbfMolecule"]
#			molecule = line["Molecule"]
			functional = line["Functional"]
			if functional == "CAM-B3LYP":
				functional = "CAMB3LYP"
			if line["no contr"] == '':
#				print 'no data for ', molecule, functional
				continue
			if molecule != args.molecule or functional != args.functional:
				continue
			excitation = line["Excitation"]
			no_contr = line["no contr"]
			for i in range(int(no_contr)):
				phia = line["occ" + str(i + 1)][:-1]
				phib = line["virt" + str(i + 1)][:-1]
				mos.append(phia)
				mos.append(phib)
		if len(mos) > 0:
			mos = np.sort(utils.remove_duplicates(np.array(mos)))
#			print(args.molecule, args.functional, mos)
			if args.inter==1:
				make_input_interactive(mos, args.molecule, args.functional, args.ngridpoints)
			else:
				make_input(mos, args.molecule, args.functional, args.ngridpoints, args.minpoints, args.maxpoints)
