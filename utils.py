import numpy as np
import re
import scipy.integrate
import sys

def cube_to_array(filename, pattern = r'[\+\-]?\d{1,}\.?\d{0,}[eE]?[\+\-]?\d{0,}'):
	f = open(filename)
	natoms = 0.
	xlim, ylim, zlim = 0., 0., 0.
	nx, ny, nz = 0, 0, 0
	dx, dy, dz = 0., 0., 0.
	vals = []
	axes=np.array(['','',''])
	stats = {}
	for i, line in enumerate(f):
		if i < 2:
			continue
		elif i==2:
			tmp = re.findall(pattern, line)
			natoms = int(tmp[0])
			xlim, ylim, zlim = float(tmp[1]), float(tmp[2]), float(tmp[3])
		elif 3 <= i and i <=5:
			tmp = re.findall(pattern, line)
			n, dV = int(tmp[0]), np.array(tmp)[1:]
			if float(dV[0]) != 0:
				nx = n
				dx = float(dV[0])
				axes[i-3] = 'x'
			if float(dV[1]) != 0:
				ny = n
				dy = float(dV[1])
				axes[i-3] = 'y'
			if float(dV[2]) != 0:
				nz = n
				dz = float(dV[2])
				axes[i-3] = 'z'
		elif 5 < i and i <= 5+abs(natoms):
			continue
		elif i==5+abs(natoms)+1 and natoms < 0:
			continue
		else:
			tmp = re.findall(pattern, line)
			vals.append(tmp)
	f.close()
	stats["natoms"] = natoms
	stats["axes"] = axes
	stats["npoints"] = np.array([nx, ny, nz])
	stats["dx"] = np.array([dx, dy, dz])
	stats["lims"] = np.array([xlim, ylim, zlim])
#	print(stats)
	# build grid
	vals = np.array([float(i) for sublist in vals for i in sublist])
	grid = np.zeros((vals.size, 3))
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = k + j*nz + i*nz*ny
				grid[x,0] = xlim + i*dx
				grid[x,1] = ylim + j*dy
				grid[x,2] = zlim + k*dz
	return grid, vals, stats

def integrate(grid, vals):
	''' integrate on a grid an array of values.
	'''
	I = scipy.integrate.simps(vals, grid)
	return I

def remove_duplicates(array):
	''' auxiliary function for the 3D integration which removes duplicates from an array
	parameters:
	array: input array
	returns:
	array: input array but with duplicates removed.
	'''
	vals = []
	vals.append(array[0])
	for i in range(1, array.shape[0]):
		found = False
		for j in range(len(vals)):
			if array[i] == vals[j]:
				found = True
		if not found:
			vals.append(array[i])
	return np.array(vals)

def integrate3D(grid, vals):
	''' integrates a function on a 3D grid
	parameters:
	grid: integration grid of the format ((x1, y1, z1), ..., (xN, yN, zN)) shape (N, 3)
	vals: values of the function to be integrated at the grid points shape (N, 1)
	returns:
	integral of the function on the input grid
	'''
	x = remove_duplicates(grid[:,0])
	y = remove_duplicates(grid[:,1])
	z = remove_duplicates(grid[:,2])
	vals = vals.reshape((x.size, y.size, z.size))
	res = np.zeros(x.size)
	dx = x[1]-x[0]
	dy = y[1]-y[0]
	dz = z[1]-z[0]
	I3 = np.sum(vals)*dx*dy*dz
	return I3

