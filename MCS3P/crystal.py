#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:32:09 2020

@author: Tobias Koehler

Crystal building module

Generate a nanocrystal from input CIF file with vacancies,cube or sphere shape
and antiphase boundaries.

Classes:
    
    Atom
    Crystal
    
"""
import numpy as np
from itertools import product

class Atom():
    """
    Class to store atom information.
    
    Attributes
    ----------
    coordinates : ndarray
        Fractional x,y,z coordinates of the atom.
    element : str
        String label of the atoms element.
    label : str
        Position label, e.g. 'Fe1'.
    probability : float
        Selection probability if a gradient is applied for vacancy generation.
    """
    __slots__ = {"coordinates", "element", "label", "probability", "isAPB", "ucn"}
    
    def __init__(self, coordinates=(0.0,0.0,0.0), element="Fe", label = None, isAPB=0, ucn=0):
        self.coordinates = np.array(coordinates)
        self.element = element
        self.label = label
        self.probability = 0
        self.isAPB = isAPB
        self.ucn = ucn
        
               
class Crystal():
    def __init__(self, settings):
        np.random.seed(settings.seed)
        
        self.diameter = settings.diameter
        self.radius = settings.diameter/2.0
        
        valid_shapes = ["Sphere", "Cube"]
        if settings.shape not in valid_shapes:
            raise ValueError("The input shape is not valid.")
        else:
            self.shape = settings.shape
        
        self.offset = settings.offset
        self.coordinates = settings.unitcell.positions
        self.elements = settings.unitcell.elements
        self.atom_numbers = {element: 0 for element in self.elements}
        self.labels = settings.unitcell.labels
        self.numAPBs = settings.numAPBs
        self.lattice_a = settings.unitcell.lattice_a
        self.lattice_b = settings.unitcell.lattice_b
        self.lattice_c = settings.unitcell.lattice_c
        self.occupancies = settings.occupancies
        self.filename = settings.filename
        self.atoms = []
        
        self.createParticle()
        
    def createParticle(self):
        self.generateAtoms()
        if self.occupancies:
            self.generateVacancies()
        if self.numAPBs:
            self.generateAPBs()
        self.generateShape()
        self.saveStructureToFile()
            
    def generateAtoms(self):
        for atomCoordinates, element,label in zip(self.coordinates, self.elements, self.labels):
            for f in product(range(self.diameter+1), repeat=3):
                self.atoms.append(Atom(coordinates= atomCoordinates + np.array(f),
                                       element = element, label = label, ucn = f))
                self.atom_numbers[element] += 1
        self.atoms = np.array(self.atoms)
            
    def gaussianDistribution(x,sig):
        n = 1/(np.sqrt(2*np.pi)*sig)
        e = np.exp(-(x/sig)**2 / 2)
        return n*e
    
    def lorentzianDistribution(x,gamma):
        return (1/np.pi)*(gamma/2)/((x**2)+(gamma/2)**2)
    
    def probabilityGradient(self,gradient_sig, plot=False):
        xp = np.arange(0,self.diameter,0.01)
        g = self.lorentzianDistribution(xp,gradient_sig)
        fp = (-1*g)+g.max()
        for atom in self.atoms:
            atom.coordinates -= self.radius
            d = np.linalg.norm(atom.coordinates)
            atom.probability = np.interp(d, xp,fp)
            atom.coordinates += self.radius
   
    def generateVacancies(self, gradientSigma=None):
        uniqueLabels = np.unique(np.array([atom.label for atom in self.atoms]))
        indices_dict = {label: [] for label in uniqueLabels}
        vacancies_dict = {label: [] for label in uniqueLabels}
                
        for index, atom in enumerate(self.atoms):
            indices_dict[atom.label].append(index)
                
        for key, indexList in indices_dict.items():
            if key in self.occupancies.keys():
                if gradientSigma:
                    self.gradient(gradient_sig = gradientSigma, plot = True)
                    p = np.array([self.atoms[index].probability for index in indexList])
                    p /= p.sum()
                    vacancies_dict[key] = np.random.choice(indexList,
                                                    int((1-self.occupancies[key])*len(indexList)),
                                                    replace = False, p=p)
                else:
                    vacancies_dict[key] = np.random.choice(indexList,
                                                    int((1-self.occupancies[key])*len(indexList)),
                                                    replace = False)
                                                    
        vacancies = np.array([vacancies_dict[key] for key in vacancies_dict.keys()])
        for vac in vacancies:
            self.atoms = np.delete(self.atoms, vac)
            
    def generateAPBs(self):
        """
        Generate antiphase-boundary.
        
        Generate an antiphase-boundary through the center of the particle. First
        set the particle into the origin, then for all atoms on one side of 
        the space diagonal, i.e. atoms whose x coordinate is larger than the y
        coordinate, get shifted along the APB by one quarter of a unit cell.
        Finally the particle is shifted back to the original position.
        """
        if self.numAPBs == 2:
            for atom in self.atoms:
                atom.coordinates -= self.radius
                if (atom.coordinates[0]- atom.coordinates[1] > (-np.sqrt(2))):
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                if (atom.coordinates[0]- atom.coordinates[1] > (np.sqrt(2))):
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                atom.coordinates += self.radius
                 
        if self.numAPBs == 1:
            for atom in self.atoms:
                atom.coordinates -= self.radius
                if (atom.coordinates[0]-atom.coordinates[1]) > 0.0:
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                if (((atom.coordinates[0]-atom.coordinates[1] < 0.0) and (atom.coordinates[0]-atom.coordinates[1] > -0.4))
                    or ((atom.coordinates[0]-atom.coordinates[1] < 0.4) and (atom.coordinates[0]-atom.coordinates[1] > 0.0))):
                    atom.isAPB = 1
                atom.coordinates += self.radius
            
    def cut_sphere(self):
        cut_crystal = []
        for atom in self.atoms:
            atom.coordinates -= self.radius
            atom.coordinates -= self.offset
            if (np.dot(atom.coordinates, atom.coordinates) < self.radius**2):
                atom.coordinates += self.radius
                atom.coordinates += self.offset
                cut_crystal.append(atom)
        self.atoms = cut_crystal
        
    def cut_cube(self, edgelength):
        cut_crystal = []
        for atom in self.atoms:
            atom.coordinates -= (edgelength/2)
            if (abs(atom.coordinates[0]) <= (edgelength/2) 
                and abs(atom.coordinates[1]) <= (edgelength/2)
                and abs(atom.coordinates[2]) <= (edgelength/2)):
                atom.coordinates += (edgelength/2)
                cut_crystal.append(atom)
        self.atoms = cut_crystal
        
    def generateShape(self):
        if self.shape == "Sphere":
            self.cut_sphere()
            volume = 4/3*np.pi*(self.diameter*self.lattice_a/10/2)**3
            self.diameter = self.diameter*self.lattice_a/10
            print("   - sphere shape generated:")
            print("\t\tdiameter: %.1f nm \n\t\tvolume: %.2f nm^3" %(self.diameter, volume))
        elif self.shape == "Cube":
            edgelength = (4/3*np.pi*(self.diameter/2)**3)**(1/3)
            self.diameter = edgelength*self.lattice_a/10
            volume = edgelength**3
            self.cut_cube(edgelength)
            print("   - cube shape generated:")
            print("\t\tedge length: %.1f unit cells, volume: %.2f" %(self.diameter, volume))
            
    def saveStructureToFile(self):
        Fe_atoms = [atom for atom in self.atoms if atom.element == "Fe"]
        xyz = [atom.coordinates for atom in Fe_atoms]
        position = [1 if (atom.label == "Fe(tet)" or atom.label == "Fe1") 
                    else 0 for atom in Fe_atoms]
        isAPB = [atom.isAPB for atom in Fe_atoms]
        unit_cell = [atom.ucn for atom in Fe_atoms]
        data = np.column_stack((xyz, position, isAPB, unit_cell))
        np.savetxt(self.filename, data, fmt='%.5f')
        
    def printStats(self, vacancies_dict, indices_dict):
        print("   - vacancies generated:")
        total_num = 0
        total_vac = 0
        for key in vacancies_dict.keys():
            total_num += len(indices_dict[key])
            total_vac += len(vacancies_dict[key])
            occ_calc = (1-(len(vacancies_dict[key])/len(indices_dict[key])))
            print("\t\t%s: %d/%d -> occ: %.2f" %(key,len(vacancies_dict[key]),len(indices_dict[key]),occ_calc))
        print("\tTotal occupancy: %.2f"%(1-(total_vac/total_num)))
        
        print("Nanoparticle generated.")
        atom_el_numbers = {}
        
        for i in self.elements:
            atom_el_numbers[i] = 0
        for i in self.atoms:
            atom_el_numbers[i.element] += 1
            
        self.atom_numbers = atom_el_numbers
            
        print("\tNumber of Fe atoms: %d\n\tNumber of O atoms: %d (not considered for MC sim)"%(atom_el_numbers["Fe"], atom_el_numbers["O"]))
        
        print("\nCrystal initialized")
        for key in self.atom_numbers.keys():
            print("\t %s: %d"%(key,self.atom_numbers[key]))
