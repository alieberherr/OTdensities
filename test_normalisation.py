import numpy as np
import utils
import argparse
import csv
import os

parser = argparse.ArgumentParser(description="Computation of Lambda descriptor")
parser.add_argument('--dir', help="working directory")

#molecules = ["DMABN","PP","acene1","acene2","acene3","acene4","acene5","betadipeptide","co","dipeptide","formaldehyde",
#             "hcl","n2","polyacetylene2","polyacetylene3","polyacetylene4","polyacetylene5","tripeptide"]
#molecules = ["betadipeptide","hcl"]
molecules = ["hcl"]
functionals = ['pbe','b3lyp','camb3lyp']

# def test_orca(args):
# 	#mos_PBE = []
# 	#mos_B3LYP = []
# 	mos_PBE = [7,8,9,10,11,13,18,25,26]
# 	mos_B3LYP=[7,8,9,10,11,12,13,14,15,17,18,21,22,23,27]
# 	mos_CAMB3LYP=[7,8,9,10,11,12,13,14,15,17,18,19,21,22,23,24,26,27,28,31,33,34]#,36]
# 
# #	orbitals = 'interactive'
# #	orbitals = 'orcainput'
# 	orbitals = 'large'
# #	orbitals = 'huge'
# 
# 	print "#-------------------------#"
# 	print "#orbitals: %s"%orbitals
# 	print "#-------------------------#"
# 
# 	print "PBE functional:"
# 	for mo in mos_PBE:
# 		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_PBE.mo%sa.cube"%(orbitals,mo))
# 		I1 = utils.integrate3D(grid,vals**2)
# 		print "mo=%s - normalisation %f"%(mo, I1)
# 	print
# 	print "B3LYP functional:"
# 	for mo in mos_B3LYP:
# 		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_B3LYP.mo%sa.cube"%(orbitals,mo))
# 		I1 = utils.integrate3D(grid,vals**2)
# 		print "mo=%s - normalisation %f"%(mo, I1)
# 	print
# 	print "CAM-B3LYP functional:"
# 	for mo in mos_CAMB3LYP:
# 		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_CAMB3LYP.mo%sa.cube"%(orbitals,mo))
# 		I1 = utils.integrate3D(grid,vals**2)
# 		print "mo=%s - normalisation %f"%(mo, I1)
        
def test_turbomole(args):
	for molecule in molecules:
		print("====",molecule,"======")
		for functional in functionals:
			print("====",functional,"======")
			if os.path.isfile(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional)):
				dict = np.load(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional),allow_pickle=True).item()
			else:
				print("failed:",molecule, functional)
				continue
			for key in dict.keys():
				phia = key.replace(" ","")
				if ("e" in phia):
					for i in ["1","2"]:
						filename = args.dir + molecule + "/" + functional + "/" + phia + i + ".cub"
						grid, vals, stats = utils.cube_to_array(filename)
						I = utils.integrate3D(grid, vals**2)
						print(f"%s%s: int |phi|^2 dV=%f"%(key,i,I))
				else:
					filename = args.dir + molecule + "/" + functional + "/" + phia + ".cub"
					grid, vals, stats = utils.cube_to_array(filename)
					I = utils.integrate3D(grid, vals**2)
					print(f"%s: int |phi|^2 dV=%f"%(key,I))
		print("==========================================")
                
                

if __name__ == '__main__':
	args = parser.parse_args()
    
	test_turbomole(args)
