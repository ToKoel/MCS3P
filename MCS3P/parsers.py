#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:45:51 2020

@author: tobiaskohler

"""        
class CifParser():
    """ 
    Class to retrieve information from input CIF files. Currently only 
    VESTA Cifs are supported. 
    
    """
    def __init__(self, filename):
        self.filename = filename
        file = open(filename, "r")
        lines = file.readlines()
        space = []
        symmetry_elements = []
        unitcell = []
        coordinatelist = []
        expanded = []
        elementlist = []
        expanded_elements = []
        labels = []
        expanded_labels = []
        self.lattice_a = 0.0
        self.lattice_b = 0.0
        self.lattice_c = 0.0
        
        for n,line in enumerate(lines):
            if "cell_length_a" in line:
                self.lattice_a = float(line.split()[1])
            elif "cell_length_b" in line:
                self.lattice_b = float(line.split()[1])
            elif "cell_length_c" in line:
                self.lattice_c = float(line.split()[1])            
            elif "space_group_name" in line:
                self.spcgr = line.split()[1]
            elif line == "\n":
                space.append(n)
                
        for line in lines[space[2]+3:space[3]]:
            l = line.strip()
            l = l[1:-1]
            symmetry_elements.append(l)
            
        for line in lines[space[3]+10:]:
            unitcell.append(line.split())
            
        for i in unitcell:
            coordinatelist.append([float(i[2]),float(i[3]),float(i[4])])
            elementlist.append(i[-1])
            labels.append(i[0])

        for k,i in enumerate(coordinatelist):
            for n in symmetry_elements:
                x = i[0]
                y = i[1]
                z = i[2]
                
                # calculate symmetry equivalent positions
                symm = n.split(',')
                new_x = eval(symm[0])
                new_y = eval(symm[1])
                new_z = eval(symm[2])
                 
                if new_x < 0:
                    new_x = 1+new_x
                if new_y < 0:
                    new_y = 1+new_y
                if new_z < 0:
                    new_z = 1+new_z
                    
                if new_x == -0.0:
                    new_x = 0.0
                if new_y == -0.0:
                    new_y = 0.0
                if new_z == -0.0:
                    new_z = 0.0
                    
                if new_x >= 1.0:
                    new_x -= 1.0
                if new_y >= 1.0:
                    new_y -= 1.0
                if new_z >= 1.0:
                    new_z -= 1.0
                    
                new_x = round(new_x,4)
                new_y = round(new_y,4)
                new_z = round(new_z,4)
                    
                expanded.append([new_x, new_y, new_z])
                expanded_elements.append(elementlist[k])
                expanded_labels.append(labels[k])
          
        
        self.positions = []
        self.elements = []
        self.labels = []
        for n,i in enumerate(expanded):
            if i not in self.positions:
                self.positions.append(i)
                self.elements.append(expanded_elements[n])
                self.labels.append(expanded_labels[n])
                
        for n,i in enumerate(self.positions):
            for k,l in enumerate(self.positions):
                if i == l and n != k:
                    print("duplicate at %d, %d"% (n,k))
                    print(i, l)
        
        print("\nCif import succesful\n")
