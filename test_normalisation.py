import numpy as np
import utils
import argparse
import csv
import os
import glob

parser = argparse.ArgumentParser(description="Test normalisation")
parser.add_argument('--dir', help="working directory")

molecules = ["DMABN","PP","acene1","acene2","acene3","acene4","acene5","betadipeptide","co","dipeptide","formaldehyde",
             "hcl","n2","polyacetylene2","polyacetylene3","polyacetylene4","polyacetylene5","tripeptide"]
molecules = ["n2","formaldehyde","co"]
molecules = ["formaldehyde"]
functionals = ['b3lyp','camb3lyp']

        
def test_turbomole(args):
	# volume integrals are written to the file wfile
	with open("../collection/normalisation.csv",'w') as wfile:
		writer = csv.DictWriter(wfile, fieldnames=["Molecule", "Functional", "Orbital", "Integral"])
		writer.writeheader()
		for molecule in molecules:
			print("====",molecule,"======")
			for functional in functionals:
				print("====",functional,"======")
				for fnam in glob.glob(args.dir+molecule+"/"+functional+"/*.cub"):
					phia = fnam.split("/")[-1][:-4]
					print("phia:",phia)
					filename = args.dir + molecule + "/" + functional + "/" + phia + ".cub"
					grid, vals, stats = utils.cube_to_array(filename)
					I = utils.integrate3D(grid, vals**2)
					print(f"%s: int |phi|^2 dV=%f"%(phia,I))
					writer.writerow({"Molecule": molecule, "Functional": functional, "Orbital": phia, "Integral": I})
					if (I < 0.99):
						print("orbital not normalised, moving on...")
						break
                

if __name__ == '__main__':
	args = parser.parse_args()
    
	test_turbomole(args)
