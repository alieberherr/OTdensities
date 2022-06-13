from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
import re

def read(filename):
    ''' reads the density from a file that is an output of ADF and returns it for plotting.
    '''
    pattern = r'[\+\-]?\d{1,}\.?\d{0,}[eE]?[\+\-]?\d{0,}'
    i = open(filename, 'r')
    i.readline()
    i.readline()
    line0 = i.readline()
    linex = i.readline()
    liney = i.readline()
    linez = i.readline()
    NA = int(re.findall(pattern, line0)[0])
    x0 = float(re.findall(pattern, line0)[1])
    y0 = float(re.findall(pattern, line0)[2])
    z0 = float(re.findall(pattern, line0)[3])
    Nx = int(re.findall(pattern, linex)[0])
    vxx = float(re.findall(pattern, linex)[1])
    Ny = int(re.findall(pattern, liney)[0])
    vyy = float(re.findall(pattern, liney)[2])
    Nz = int(re.findall(pattern, linez)[0])
    vzz = float(re.findall(pattern, linez)[3])
    for l in range(abs(NA)):
        i.readline()
    if NA < 0:
        i.readline()
    alldata = []
    for line in i:
        dat = re.findall(pattern, line)
        for entry in dat:
            alldata.append(float(entry))
    alldata = np.array(alldata)
    alldata = alldata.reshape((Nz, Ny, Nx))
    i.close()
    return alldata, Nx, Ny, Nz

def plot_xy(density, z0):
    nx, ny  =density.shape[2], density.shape[1]
    x, y = np.arange(nx) - 0.5*nx, np.arange(ny) - 0.5*ny
    X, Y = np.meshgrid(x, y)
    D = density[z0,:,:]
    #plot(X = X, Y= Y,Z =  D**2)
    plt.contourf(X,Y,D**2)
    #plt.show()
    plt.savefig(fileout + "_xy.png")

def plot_xz(density, y0):
    nx, nz = density.shape[2], density.shape[0]
    x, z = np.arange(nx) - 0.5*nx, np.arange(nz) - 0.5*nz
    X, Z = np.meshgrid(x, z)
    D = density[:,y0,:]
    #plot(X = X, Y = Z, Z =  D**2)
    plt.contourf(X,Z,D**2)
    #plt.show()
    plt.savefig(fileout + "_xz.png")

def plot_yz(density, x0):
    ny, nz = density.shape[1], density.shape[0]
    y, z  =np.arange(ny) - 0.5*ny, np.arange(nz) - 0.5*nz
    Y, Z = np.meshgrid(y, z)
    D = density[:,:,x0]
    #plot(X = Y, Y = Z, Z = D**2)
    plt.contourf(Y,Z,D**2)
    #plt.show()
    plt.savefig(fileout + "_yz.png")

def plot(X, Y, Z):
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.plot_surface(X, Y, Z)

if __name__ == '__main__':
    filename = sys.argv[1]
    fileout = sys.argv[2]
    density, Nx, Ny, Nz = read(filename)
    plot_xy(density, int(Nz/2))
    plt.clf()
    plot_yz(density, int(Nx/2))
    plt.clf()
    plot_xz(density, int(Ny/2))
    plt.clf()
