import numpy as np
import utils


if __name__ == '__main__':
	mos_PBE = [7,8,9,10,11,13,18,25,26]
	mos_B3LYP=[7,8,9,10,11,12,13,14,15,17,18,21,22,23,27]
	mos_CAMB3LYP=[7,8,9,10,11,12,13,14,15,17,18,19,21,22,23,24,26,27,28,31,33,34]#,36]

#	orbitals = 'interactive'
#	orbitals = 'orcainput'
	orbitals = 'large'
#	orbitals = 'huge'

	print "#-------------------------#"
	print "#orbitals: %s"%orbitals
	print "#-------------------------#"

	print "PBE functional:"
	for mo in mos_PBE:
		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_PBE.mo%sa.cube"%(orbitals,mo))
		I1 = utils.integrate3D(grid,np.abs(vals))
		print "mo=%s - absolute value %f"%(mo, I1)
	print
	print "B3LYP functional:"
	for mo in mos_B3LYP:
		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_B3LYP.mo%sa.cube"%(orbitals,mo))
		I1 = utils.integrate3D(grid,np.abs(vals))
		print "mo=%s - absolute value %f"%(mo, I1)
	print
	print "CAM-B3LYP functional:"
	for mo in mos_CAMB3LYP:
		grid, vals = utils.cube_to_array("/cluster/scratch/anninal/OTproject/case_normalisation/formaldehyde/%s/formaldehyde_CAMB3LYP.mo%sa.cube"%(orbitals,mo))
		I1 = utils.integrate3D(grid,np.abs(vals))
		print "mo=%s - absolute value %f"%(mo, I1)
