import numpy as np
import utils
import argparse
import csv
import os

# calculate < r^2 > for all orbitals in the data set

def rsquared(phi):
	grid, vals, stats = utils.cube_to_array(phi)
	xc = .5*(grid[-1,0] + grid[0,0])
	yc = .5*(grid[-1,1] + grid[0,1])
	zc = .5*(grid[-1,2] + grid[0,2])

	x = utils.remove_duplicates(grid[:,0])
	y = utils.remove_duplicates(grid[:,1])
	z = utils.remove_duplicates(grid[:,2])
	dx = x[1]-x[0]
	dy = y[1]-y[0]
	dz = z[1]-z[0]
	rc = np.sqrt(xc**2 + yc**2 + zc**2)
	
	rgrid = np.sqrt(np.sum(grid**2, axis=1))
	I = np.sum(vals**2 * (rgrid - rc)**2) * dx * dy * dz
	return I

if __name__=='__main__':
	pass
#	phi = "/u/dem/chem1614/Documents/projects/optimal-transport-excitations/molecules_turbomole/acene1/pbe/1au.cub"
#	rsquared(phi)
	dir = "/u/dem/chem1614/Documents/projects/optimal-transport-excitations/collection/calc_turbomole"

	f = open("../collection/rsquared/rsquared.csv",'w')
	f.write("Molecule,Functional,Orbital,<r^2>\n")
	molecules = os.listdir(dir)
	print(molecules)
	for molecule in molecules:
		if os.path.isdir("%s/%s"%(dir,molecule)):
			if molecule == "basissets":
				continue
			print(molecule)
			for functional in ['pbe','b3lyp','camb3lyp']:
				for phi in os.listdir("%s/%s/%s/orbs_8/"%(dir,molecule,functional)):
					if not phi.endswith(".cub"):
						continue
					fnam = dir + "/" + molecule + "/" + functional + "/orbs_8/" + phi
					val = rsquared(fnam)
					f.write("%s,%s,%s,%f\n"%(molecule,functional,phi,val))
