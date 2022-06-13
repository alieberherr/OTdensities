import ot
import numpy as np
import time
import utils

n=100
eps=0.5
x=np.arange(n+1,dtype=np.float64)[1:]
a = np.arange(1,n+1)/np.sum(x)
b = np.arange(1,n+1)[::-1]/np.sum(x)

# filenamea = "test/polyacetylene2/pbe/orbs_adjgrid_spacing0_8/2au.cub"
# filenameb = "test/polyacetylene2/pbe/orbs_adjgrid_spacing0_8/1bg.cub"
# 
# grid,vals1,stats1 = utils.cube_to_array(filenamea)
# grid,vals2,stats2 = utils.cube_to_array(filenameb)
# 
# vals1 = vals1/np.sqrt(np.sum(vals1**2))
# vals2 = vals2/np.sqrt(np.sum(vals2**2))

print("Unstabilised POT Sinkhorn:")
start = time.time()
M = ot.dist(x.reshape((n,1)),x.reshape((n,1)))
# M = ot.dist(grid,grid)
plan=ot.bregman.sinkhorn(a,b,M,eps,verbose=True)
# plan=ot.bregman.sinkhorn(vals1**2+1e-16,vals2**2+1e-16,M,eps,verbose=True)
print(f"Final distance: %f"%(np.sum(np.multiply(M,plan))))
end = time.time()
print(f"time: %f"%(end-start))

print("Log-stabilised POT Sinkhorn:")
start = time.time()
M = ot.dist(x.reshape((n,1)),x.reshape((n,1)))
# M = ot.dist(grid,grid)
plan=ot.bregman.sinkhorn(a,b,M,eps,verbose=True,method="sinkhorn_stabilized",numItermax=50000)
# plan=ot.bregman.sinkhorn(vals1**2+1e-16,vals2**2+1e-16,M,eps,verbose=True)
print(f"Final distance: %f"%(np.sum(np.multiply(M,plan))))
end = time.time()
print(f"time: %f"%(end-start))
