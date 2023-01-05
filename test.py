import ot
import numpy as np
import time
import utils

n=100
eps=0.5
x=np.arange(n+1,dtype=np.float64)[1:]
a = np.arange(1,n+1)/np.sum(x)
b = np.arange(1,n+1)[::-1]/np.sum(x)

filenamea = "acene1/1au.cub"
filenameb = "acene1/2b3g.cub"
# 
grid,vals1,stats1 = utils.cube_to_array(filenamea)
grid,vals2,stats2 = utils.cube_to_array(filenameb)
# 
vals1 = vals1/np.sqrt(np.sum(vals1**2))
vals2 = vals2/np.sqrt(np.sum(vals2**2))

print("Unstabilised POT Sinkhorn:")
start = time.time()
#M = ot.dist(x.reshape((n,1)),x.reshape((n,1)))
M = ot.dist(grid,grid)
#plan=ot.bregman.sinkhorn(a,b,M,eps,verbose=True)
plan=ot.bregman.sinkhorn(vals1**2+1e-16,vals2**2+1e-16,M,eps,verbose=True)
print(f"Final distance: %f"%(np.sum(np.multiply(M,plan))))
end = time.time()
print(f"time: %f"%(end-start))

#print("Log-stabilised POT Sinkhorn:")
#start = time.time()
#M = ot.dist(x.reshape((n,1)),x.reshape((n,1)))
#M = ot.dist(grid,grid)
#plan=ot.bregman.sinkhorn(a,b,M,eps,verbose=True,method="sinkhorn_stabilized",numItermax=250000)
#plana=ot.bregman.sinkhorn(a,a,M,eps,verbose=True,method="sinkhorn_stabilized",numItermax=250000)
#planb=ot.bregman.sinkhorn(b,b,M,eps,verbose=True,method="sinkhorn_stabilized",numItermax=250000)
# plan=ot.bregman.sinkhorn(vals1**2+1e-16,vals2**2+1e-16,M,eps,verbose=True)
#print(f"Final distance: %f"%(np.sum(np.multiply(M,plan))))

#entr = 0.
#entra = 0.
#entrb = 0.
#for i in range(vals1.size):
#	for j in range(vals2.size):
#		if plan[i,j]>1e-8:
#			entr += plan[i,j]*np.log(plan[i,j])
#		if plana[i,j]>1e-8:
#			entra += plana[i,j]*np.log(plana[i,j])
#		if planb[i,j]>1e-8:
#			entrb += planb[i,j]*np.log(planb[i,j])
#print("entr. part:",eps*entr)
#print("Symmetrised cost:",np.sum(np.multiply(M,plan))+eps*entr-.5*(np.sum(np.multiply(M,plana))+eps*entra + np.sum(np.multiply(M,planb))+eps*entrb))

#end = time.time()
#print(f"time: %f"%(end-start))
