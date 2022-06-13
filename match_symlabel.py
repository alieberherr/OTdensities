import numpy as np
import re
import csv
import os

pattern_index = re.compile("\d{1,}")
pattern_sym = re.compile("\d{1,}\s[abe]\d?[gu]?[\'\"]?")

molecules = ["DMABN","PP","acene1","acene2","acene3","acene4","acene5","betadipeptide","co","dipeptide","formaldehyde",
             "hcl","n2","polyacetylene2","polyacetylene3","polyacetylene4","polyacetylene5","tripeptide"]
functionals = ['pbe','b3lyp','camb3lyp']

def make_dictionary(ex_filename):
    
    with open(ex_filename, 'r') as rfile:
        reader = csv.DictReader(rfile)
        for i,line in enumerate(reader):
            if i==0:
                continue
               
            molecule = line["\ufeffMolecule"]
            functional = line["Functional"]
            # get which 
            ncntrib = int(line["no contr"])
            if ncntrib == 0:
                continue
                
            wanted = []
            mapped = []
            for i in range(ncntrib):
                phia = line["occ" + str(i + 1)]
                phib = line["virt" + str(i + 1)]
                if not phia in wanted:
                    wanted.append(phia)
                if not phib in wanted:
                    wanted.append(phib)
                    
            print("wanted:",wanted)
                
            if os.path.isfile(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional)):
                dict = np.load(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional),allow_pickle=True).item()
            else:
                dict = {}
                
            eiger = open(f"../molecules_turbomole/%s/%s/eiger.out"%(molecule,functional),'r')
            for bline in eiger:
                tmp1 = re.findall(pattern_index,bline)
                tmp2 = re.findall(pattern_sym,bline)
                if len(tmp1) > 0 and len(tmp2) > 0 and tmp2[0] in wanted and tmp2[0] not in dict.keys():
                    dict[tmp2[0]] = tmp1[0]
                    mapped.append(tmp1[0])
            
            print(molecule, functional, ": ", dict)
            print()

            np.save(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional),dict)

def make_inputs():
    for molecule in molecules:
        for functional in functionals:
            if not os.path.isfile(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional)):
                continue
            dict = np.load(f"../data/excitations/turbomole/maps/%s_%s_map.npy"%(molecule,functional),allow_pickle=True).item()
            
            # write the wanted numbers to file
            o = open(f"../data/excitations/turbomole/maps/%s_%s_map.txt"%(molecule,functional),'w+')
            for curr in dict.keys():
                o.write(f"%s "%dict[curr])
            o.close()

make_dictionary("../data/excitations/turbomole/ex_all_turbomole.csv")
make_inputs()
