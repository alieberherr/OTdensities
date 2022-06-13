#!/bin/bash


bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part1.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part1.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part2.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part2.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part3.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part3.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part4.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part4.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part5.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part5.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part6.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part6.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part7.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part7.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part8.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part8.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part9.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part9.csv --orbitals orbs_adjgrid_spacing0_8 " 

bsub -W 100:00 -R "rusage[mem=400000]" "julia theta_new.jl --resultsfile data/theta/div_10cores/res_only_daug_part10.csv --exfile data/excitations/turbomole/div_10cores/ex_only_daug_part10.csv --orbitals orbs_adjgrid_spacing0_8 " 

