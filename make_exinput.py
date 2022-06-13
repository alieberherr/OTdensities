import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--ncores",type=int,help="number of cores wanted")


if __name__=="__main__":
    args = parser.parse_args()
    mode = "without_daug"
    # mode = "only_daug"
    
    all_ex = open(f"../data/transitions/turbomole/ex_%s.csv"%(mode),"r")
    num_lines = sum(1 for line in all_ex)-1
    all_ex.close()
    all_ex = open(f"../data/transitions/turbomole/ex_%s.csv"%(mode),"r")
    header = all_ex.readline()
    print(header)
    
    blocksize=num_lines//args.ncores
    finalblock=num_lines-args.ncores*blocksize
    blocks=np.zeros(args.ncores,dtype=int)
    for i in range(args.ncores-1):
        blocks[i] = blocksize
    
    if finalblock == 0:
        print(args.ncores," times ","blocksize=",blocksize,"\n")
        blocks[-1] = blocksize
    else:
        print(args.ncores-1," times ","blocksize=",blocksize,"\n")
        print("finalblock=",finalblock+blocksize,"\n")
        blocks[-1] = finalblock+blocksize
    
    print("blocks=",blocks)
    for i in range(args.ncores):
        ex_curr = open(f"../data/transitions/turbomole/div_%icores/ex_%s_part%i.csv"%(args.ncores,mode,i+1),"w")
        ex_curr.write(header)
        for j in range(blocks[i]):
            line = all_ex.readline()
            ex_curr.write(line)
    
        ex_curr.close()
    all_ex.close()

    driver = open(f"driver_%s_%icores.sh"%(mode,args.ncores),"w")
    driver.write("#!/bin/bash\n\n\n")
    
    for i in range(args.ncores):
        driver.write(f"bsub -W 100:00 -R \"rusage[mem=400000]\" \"julia theta_new.jl --resultsfile data/theta/div_%icores/res_%s_part%i.csv --exfile data/excitations/turbomole/div_%icores/ex_%s_part%i.csv --orbitals orbs_adjgrid_spacing0_8 \" \n\n"%(args.ncores,mode,i+1,args.ncores,mode,i+1))
    driver.close()
    
    
    
