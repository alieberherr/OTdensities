import numpy as np
import utils

filename = "../molecules_turbomole/co/pbe_test/7a1.cub"
grid, vals, stats = utils.cube_to_array(filename)

I = utils.integrate3D(grid, vals**2)
print(filename,":")
print(f"int |phi|^2 dV=%f"%(I))

# filename = "../molecules_turbomole/co/pbe/6a1.cub"
# grid, vals, stats = utils.cube_to_array(filename)
# 
# I = utils.integrate3D(grid, vals**2)
# print(filename,":")
# print(f"int |phi|^2 dV=%f"%(I))
