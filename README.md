# OTdensities

Codes and data files to accompany the report "Optimal Transport distances to characterise electronic excitations", currently found at arXiv:2308.09118.

# Code

| name | description |
|------|-------------|
| lambda.py | calculates $\Lambda$ diagnostic | 
| lambda_sq.py | calculates modified $\Lambda$ (same as $S$ in eq. [...] from [2]) |
| normalise_theta.py | calculate $\Theta'$ from pre-computed $\Theta$ and <r<sup>2</sup>> |
| rsquared.py | calculate  <r<sup>2</sup>> for a set of orbitals |
| theta.py | calculates $\Theta$ diagnostic |
| theta_conv.py | calculates $\Theta$ diagnostic for decreasing regularisation parameters |
| utils.py | general functions such as reading a .cub file into an array |

<!--# Data

1. Input files: Input files for the TDDFT calculations with Turbomole and Orca
2. Excitation energies: Values of TDDFT excitation energies from Turbomole, Orca and article [1]
3. $\Lambda$ (absolute value of orbitals as well as absolute value squared) and $\Theta$ values

-->


-----------------------
[1] Peach et al., *J. Chem. Phys*, **128**, (2008)

[2] Giesbertz et al. (???)
