#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:59:36 2020

@author: tobiaskohler
"""

import crystal
import parsers
import numpy as np
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=str)
    parser.add_argument('--occ_oct', type=float)
    parser.add_argument('--occ_tet', type=float)
    parser.add_argument('--diameter', type=int)
    parser.add_argument('--number', type=int)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--APB', action="store_true")
    parser.add_argument('--template', action="store_true")
    parser.add_argument('--ciffile', type=str)
    parser.add_argument('--atomlabels', type=str)
    parser.add_argument('--occupancies', type=str)
    parser.add_argument('--shape', type=str)
    parser.add_argument('--latticepar', type=str)
    
    args = parser.parse_args()
    
    #Fe3O4_Fd3m = "/Users/tobiaskohler/PhD/Debye-TK/cif-files/Fd-3m_perfect.cif" 
    #cif_file = Fe3O4_Fd3m
    
    cif_file = args.ciffile
    shape = "Sphere"
    shape = args.shape
    n_APBs = 1
    offset = np.array([0.5,0.5,0.5])
    gradient = None
    
    APB = args.APB
    template = args.template
    path = args.output
    occ_oct = args.occ_oct
    occ_tet = args.occ_tet
    diameter = args.diameter
    number = args.number
    seed = args.seed
    
    np.random.seed(seed)
    rands = np.random.randint(low=16430104, high=20210503, size=100)
    
    labels = (args.atomlabels).split(",")
    occs = (args.occupancies).split(",")
    occupancies = {}
    for l,o in zip(labels,occs):
        occupancies[l] = float(o)
        
    latticepars = (args.latticepar).split(",")
        
    #occupancies = {'Fe(oct)':occ_oct, 'Fe(tet)':occ_tet}
    

    if occupancies == {}:
        occ = False
    else:
        occ = True

    saveStructure = path + "D%d_structure_" %(diameter)

    if APB:
        saveStructure += "APB"
    else:
        saveStructure += "noAPB"


    if template:
        saveStructure += "_template"
        SEED = rands[0]
    else:
        saveStructure += "_%d"%(number)
        SEED = rands[number]

    unitcell = parsers.CifParser(cif_file)
    unitcell.lattice_a = float(latticepars[0])
    unitcell.lattice_b = float(latticepars[1])
    unitcell.lattice_c = float(latticepars[2])
        
    Crystal = crystal.Crystal(diameter, unitcell, occupancies, shape, n_APBs, offset)
        
    Crystal.build_nanoparticle(APB = APB,
                           occ = occ,
                           gradient = gradient,
                           plot = False,
                           SEED = SEED)
        
    Crystal.output_crystal_structure(saveStructure, oxygen=False)

    print("----------------------------\n")
    
if __name__ == "__main__":
    main()
