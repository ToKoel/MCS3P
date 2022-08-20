#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:59:36 2020

@author: tobiaskohler
"""

from crystal import Crystal
import parsers
import numpy as np
import argparse

class Settings:
    def __init__(self, diameter, shape, offset, fileName, occupancies, seed, numAPBs, unitcell):
        self.diameter = diameter
        self.shape = shape
        self.offset = offset
        self.filename = fileName
        self.occupancies = occupancies
        self.seed = seed
        self.numAPBs = numAPBs
        self.unitcell = unitcell
        
def generateFileName(args):
    fileName = args.output + f"/D{args.diameter}_structure_"
    if args.numAPBs > 0:
        fileName += "APB"
    else:
        fileName += "noAPB"

    if args.template:
        fileName += "_template"
    else:
        fileName += f"_{args.number}"
    return fileName

def generateParticle(args):
    np.random.seed(args.seed)
    rands = np.random.randint(low=16430104, high=20210503, size=100)
    if args.template:
        SEED = rands[0]
    else:
        SEED = rands[args.number]
    
    labels = (args.atomlabels).split(",")
    occupancies = (args.occupancies).split(",")
    occupancyDict = {label: float(occupancy) for (label, occupancy) in zip(labels, occupancies)}

    unitcell = parsers.CifParser(args.ciffile)
    latticepars = (args.latticepar).split(",")
    unitcell.lattice_a = float(latticepars[0])
    unitcell.lattice_b = float(latticepars[1])
    unitcell.lattice_c = float(latticepars[2])
    
    offset = np.array([0.5,0.5,0.5])
    
    Crystal(Settings(args.diameter,
                     args.shape,
                     offset,
                     generateFileName(args),
                     occupancyDict,
                     SEED,
                     args.numAPBs, unitcell))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=str)
    parser.add_argument('--occ_oct', type=float)
    parser.add_argument('--occ_tet', type=float)
    parser.add_argument('--diameter', type=int)
    parser.add_argument('--number', type=int)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--numAPBs', type=int, default=0)
    parser.add_argument('--template', action="store_true")
    parser.add_argument('--ciffile', type=str)
    parser.add_argument('--atomlabels', type=str)
    parser.add_argument('--occupancies', type=str)
    parser.add_argument('--shape', type=str)
    parser.add_argument('--latticepar', type=str)
    args = parser.parse_args()
    
    generateParticle(args)
