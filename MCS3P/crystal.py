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
import sys


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
    
    def __init__(self, coordinates=(0.0,0.0,0.0), element="Fe", label = None,isAPB=0, ucn=0):
        self.coordinates = np.array(coordinates)
        self.element = element
        self.label = label
        self.probability = 0
        self.isAPB = isAPB
        self.ucn = ucn
        
               
class Crystal():
    """
    Class to store crystal information and perform calculations to build the 
    nanocrystal.
    
    Attributes
    ----------
    diameter : float
        The nanocrystal diameter, in case of cubic shape the edge length.
    shape : str
        Shape of the nanocrystal ("Cube" or "Sphere").
    radius : float
        Half of the diameter.
    atoms : list
        List of Atom objects.
    unitcell : ndarray
        Array of atomic coordinates in the unitcell.
    elements : chararray
        Array of element labels for atoms in the unitcell.
    labels : chararray
        Array of site labels for atoms in the unitcell.
    cif_filename : str
        Filename of the input CIF file.
    atom_numbers : dictionary
        Number of atoms per element.
    lattice_a : float
        Lattice parameter a of unitcell in Angstrom.
    lattice_b : float
        Lattice parameter b of unitcell in Angstrom.
    lattice_c : float
        Lattice parameter c of unitcell in Angstrom.
    occupancies : dictionary
        Occupancies for each symmetrically different iron position.
    
    Methods
    -------
    gaussian(x,sig) 
        Generate gaussian distribution around 0.
    lorentzian(self,x,gamma)   
        Generate lorentzian distribution around 0.
    gradient(gradient_sig, plot=False)
        Calculate selection probabilities for atoms based on a gradient.
    generate_vacancies(occ,gradient=None, SEED=0)
        Introduce vacancies on iron positions in the nanocrystal.
    generate_APB()
        Generate 1/4[110]-APB. 
    cut_sphere(diameter)
        Cut the nanocrystal into spherical shape.
    cut_cube(edgelength)
        Cut the nanocrystal into cube shape.
    return_coordinate_list()
        Return the list of all atom coordinates.
    calculate_distance_array()
        Calculate the pair distances in the nanocrystal.
    calculate_formfactor_matrix()
        Setup the matrix containing string lables for formfactor products.
    build_nanoparticle(APB, occ,gradient, plot)
        Build the nanoparticle with different options.
    get_properties()
        Get the distance, f_ij and f_ii arrays.
    rotate(alpha,beta,gamma)
        Rotate the crystal.
        
    """
    def __init__(self, diameter, unitcell,occupancies, shape, n_APBs, offset):
        """
        Initialize a crystal object. 
        
        Called automatically on initialization of a crystal object. The unit
        cell contents are provided as parameters, and the crystal is built up 
        by repeating the unit cell in all dimensions n times, where n is 
        specified by the parameter diameter. In the repitiion loop Atom objects
        are initialized and stored in self.Atoms. self.Atoms then contains a 
        list of all Atoms in the nanocrystal.

        Parameters
        ----------
        diameter : Int
            Number of unit cell repititions in x,y and z.
        unitcell : unitcell object
            Object containing the unitcell information obtained from CifParser.
        shape : str
            If "Sphere" cut the crystal in spherical shape with diameter 
            specified in crystal initialization.
            If "Cube" cut into cube shape with same volume as a sphere of diameter 
            specified in crystal initialization.
       

        Returns
        -------
        None.

        """
        self.diameter = diameter
        valid_shapes = ["Sphere", "Cube"]
        if shape not in valid_shapes:
            raise ValueError("The input shape is not valid.")
        else:
            self.shape = shape
        self.radius = diameter/2.0
        self.offset = offset
        self.atoms = list()
        self.unitcell = unitcell.positions
        self.elements = unitcell.elements
        self.labels = unitcell.labels
        self.n_APBs = n_APBs
        self.cif_filename = ((unitcell.filename).split('/'))[-1]
        self.atom_numbers = {}
        self.lattice_a = unitcell.lattice_a
        self.lattice_b = unitcell.lattice_b
        self.lattice_c = unitcell.lattice_c
        self.occupancies = occupancies
        self.formfactor_array = []
        self.a_numbers = 0
        
        
        for i in self.elements:
            self.atom_numbers[i] = 0
        
        for a,e,l in zip(self.unitcell, self.elements, self.labels):
            for f in product(range(self.diameter+1), range(self.diameter+1), range(self.diameter+1)):
                self.atoms.append(Atom(coordinates= a + np.array(f), 
                                            element = e, label = l, ucn = f))
                self.atom_numbers[e] += 1
                
        self.atoms = np.array(self.atoms)
        print("\nCrystal initialized")
        for key in self.atom_numbers.keys():
            print("\t %s: %d"%(key,self.atom_numbers[key]))
            
    def gaussian(self,x,sig): 
        """ Gaussian distribution """
        n = 1/(np.sqrt(2*np.pi)*sig)
        e = np.exp(-(x/sig)**2 / 2)
        return n*e
    
    def lorentzian(self,x,gamma):
        """ Lorentzian distribution """
        return (1/np.pi)*(gamma/2)/((x**2)+(gamma/2)**2)
    
    def gradient(self,gradient_sig, plot=False):
        """
        Setup a probability gradient for vacancy formation.
        
        The Gradient is used to calculate a selection probability for every 
        atom in the structure.

        Parameters
        ----------
        gradient_sig : float
            Sigma parameter for distribution functions.
        plot : Bool, optional
            Plot the probabilities against the distance from particle center. 
            The default is False.

        Returns
        -------
        None.

        """
        ps = []
        ds = []
        xp = np.arange(0,self.diameter,0.01)
        #g = self.gaussian(xp,gradient_sig)
        g = self.lorentzian(xp,gradient_sig)
        fp = (-1*g)+g.max()
        #fp = g
        for atom in self.atoms:
            atom.coordinates -= self.radius
            d = np.linalg.norm(atom.coordinates)        
            ds.append(d)
            p =np.interp(d, xp,fp)
            ps.append(p)
            atom.probability = p
            atom.coordinates += self.radius
   
        
        
    def generate_vacancies(self, occ,gradient=None, SEED=0):
        """
        Generate vacancies on iron sites.
        
        This is only intended for maghemite or magnetite structures. 
        
        First the list of atom labels is generated from the list of atoms. Then
        the unique labels and their occurrances are determined with np.unique().
        Dictionaries are set up for the atom labels, the indices and the vacancies.
        After that the keys relating to oxygen in the structure are stored in 
        O_keys. Now in an iteration over all atoms the atom indices corresponding
        to the labels determined previously are stored in indices_dict. The keys
        relating to oxygen are removed. A check is performed if the input keys
        match the ones determined from the crystal. If not the program is 
        terminated. Random indices are drawn from the indices lists corresponding
        to the iron positions. The number of atoms to be removed is determined
        from the occupancy factor given in the input dictionary. The retrieved 
        indices are stored all together in vacancies_merged, which is finally
        used to remove the selected atoms from the crystal.

        Parameters
        ----------
        occ : Dictionary
            Dictionary containing the site labels with occupancies.
        gradient : Bool
            Use a gradient to determine the probability of selection.
        SEED : Int
            Seed value for numpy random.

        Returns
        -------
        Int.

        """
        np.random.seed(SEED)
        
        if gradient != None:
            self.gradient(gradient_sig=gradient,plot = True)
        
        
        label_list = []
        label_dict = {}
        indices_dict = {}
        vacancies_dict = {}
        
        
        for n,atom in enumerate(self.atoms):
            label_list.append(atom.label)
            
        label_list = np.array(label_list)
        uniq, counts = np.unique(label_list, return_counts=True)
            
        for n,l in enumerate(uniq):
            label_dict[l] = counts[n]
            indices_dict[l] = []
            vacancies_dict[l] = []
         
        O_keys = []
        for i in indices_dict.keys():
            if "O" in i:
                O_keys.append(i)
                 
        for c,a in enumerate(self.atoms):
            for n in range(len(uniq)):
                if a.label == uniq[n]:
                    indices_dict[a.label].append(c)
                    
        for i in O_keys:
            indices_dict.pop(i)
            vacancies_dict.pop(i)
                    
        for i in occ.keys():
            if i not in vacancies_dict.keys():
                print("Occupancy dict contains wrong labels!")
                print("Valid keys are: ", vacancies_dict.keys())
                print("Vacancies were not set up.")
                return 1
        
        for i in indices_dict.keys():
            if gradient:
                p = []
                for a in indices_dict[i]:
                    p.append(self.atoms[a].probability)
                p=np.array(p)
                p/=p.sum()
                vacancies_dict[i] = np.random.choice(indices_dict[i], 
                                                 int((1-occ[i])*len(indices_dict[i])), 
                                                 replace = False, p=p)
            else:
                vacancies_dict[i] = np.random.choice(indices_dict[i], 
                                                 int((1-occ[i])*len(indices_dict[i])), 
                                                 replace = False)
        vacancies = []
        vacancies_merged = []
        print("   - vacancies generated:")
        total_num = 0
        total_vac = 0
        for i in vacancies_dict.keys():
            total_num += len(indices_dict[i])
            total_vac += len(vacancies_dict[i])
            vacancies.append(vacancies_dict[i])
            occ_calc = (1-(len(vacancies_dict[i])/len(indices_dict[i])))
            print("\t\t%s: %d/%d -> occ: %.2f" %(i,len(vacancies_dict[i]),len(indices_dict[i]),occ_calc))
        print("\tTotal occupancy: %.2f"%(1-(total_vac/total_num)))
        
        for i in vacancies:
            vacancies_merged += i.tolist()
        
        self.atoms = np.delete(self.atoms, vacancies_merged)

        return 0
    
        

    def generate_APB(self,APB, n):
        """
        Generate antiphase-boundary.
        
        Generate an antiphase-boundary through the center of the particle. First
        set the particle into the origin, then for all atoms on one side of 
        the space diagonal, i.e. atoms whose x coordinate is larger than the y
        coordinate, get shifted along the APB by one quarter of a unit cell.
        Finally the particle is shifted back to the original position.
        
        Returns
        -------
        None.
        """
        if n == 2:
            for atom in self.atoms:
                atom.coordinates -= self.radius
            #    atom.coordinates -= self.offset
                if (atom.coordinates[0]- atom.coordinates[1] > (-np.sqrt(2))):
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                if (atom.coordinates[0]- atom.coordinates[1] > (np.sqrt(2))):
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                atom.coordinates += self.radius
              #  atom.coordinates += self.offset
                 
        if n == 1:
            for atom in self.atoms:
                atom.coordinates -= self.radius
                #atom.coordinates -= self.offset
                #if (atom.coordinates[0]-atom.coordinates[1]) > 0.0:
                
                if (atom.coordinates[0]-atom.coordinates[1]) > 0.0:
                    atom.coordinates += np.array([0.25, 0.25, 0.0])
                if (((atom.coordinates[0]-atom.coordinates[1] < 0.0) and (atom.coordinates[0]-atom.coordinates[1] > -0.4))
                    or ((atom.coordinates[0]-atom.coordinates[1] < 0.4) and (atom.coordinates[0]-atom.coordinates[1] > 0.0))):
                    atom.isAPB = 1
                atom.coordinates += self.radius
               # atom.coordinates += self.offset
        
            
    def cut_sphere(self, diameter,offset):
        """
        Generate sphere shape. 
        
        Cut the particle into a spherical shape by appending only atoms with
        x^2 + y^2 + z^2 < r^2 to the new structure.
        
        Parameters
        ----------
        diameter : Float
            Diameter of the nanocrystal after shaping in fractions of unit cells.

        Returns
        -------
        None.
        
        """
        radius = diameter/2.0
        cut_crystal = []
        for atom in self.atoms:
            atom.coordinates -= radius
            atom.coordinates -= offset
            if (np.dot(atom.coordinates, atom.coordinates) < radius**2):
                atom.coordinates += radius
                atom.coordinates += offset
                cut_crystal.append(atom)
        self.atoms = cut_crystal
        
    def cut_cube(self, edgelength):
        """
        Generate cube shape.
        
        Remove atoms from initialized crystal to generate a cube shaped crystal
        with specified edgelength.

        Parameters
        ----------
        edgelength : Float
            Edge length of the nanocrystal after shaping in fractions of unit 
            cells.

        Returns
        -------
        None.

        """
        cut_crystal = []
        for atom in self.atoms:
            atom.coordinates -= (edgelength/2)
            if (abs(atom.coordinates[0]) <= (edgelength/2) 
                and abs(atom.coordinates[1]) <= (edgelength/2)
                and abs(atom.coordinates[2]) <= (edgelength/2)):
                atom.coordinates += (edgelength/2)
                cut_crystal.append(atom)
        self.atoms = cut_crystal
            
  
        
    def build_nanoparticle(self, APB, occ,gradient, plot,SEED):
        """
        Setup the nanoparticle shape and calculate all parameters.

        Parameters
        ----------
        APB : Bool
            If True generate an APB along [110] through the particle center.
        occ : dictionary
            Set the occupancies iron positions. 
        gradient : Float
            Set sigma parameter of gaussian distribution for probabilities.
        plot : Bool
            If True the crystal is plotted in 3D

        Returns
        -------
        None.

        """
        print("\nBuilding Nanoparticle...")
        if occ:
            res = self.generate_vacancies(self.occupancies,gradient, SEED=SEED)
            if res == 1:
                sys.exit("terminating program")
           # elif res == 0:
           #     print("   - vacancies generated")
        
        if APB:
            self.generate_APB(APB,self.n_APBs)
            print("   - APB generated")
            
        if self.shape == "Sphere":
            self.cut_sphere(self.diameter,self.offset)
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
        
        
        print("Nanoparticle generated.")
        atom_el_numbers = {}
        
        for i in self.elements:
            atom_el_numbers[i] = 0
        for i in self.atoms:
            self.a_numbers += 1
            atom_el_numbers[i.element] += 1
            
        self.atom_numbers = atom_el_numbers
            
        print("\tNumber of Fe atoms: %d\n\tNumber of O atoms: %d (not considered for MC sim)"%(atom_el_numbers["Fe"], atom_el_numbers["O"]))
        
        
            
        
    def get_properties(self):
        """
        Return the calculated crystal properties for further use.

        Returns
        -------
        numpy.ndarray
            Array containing pair wise distances.
        numpy.ndarray
            Array containing pair wise form factor string-placeholders.
        numpy.chararray
            Array containing pair wise form factor string-placeholders of the
            diagonal elements, i.e. the entries with equal indices.

        """
        return self.distance_array, self.formfactor_array, self.formfactors_same
    
    def output_crystal_structure(self, filename, oxygen):
        Fe_atoms = [atom for atom in self.atoms if atom.element == "Fe"]
        
        x = [atom.coordinates[0] for atom in Fe_atoms]
        y = [atom.coordinates[1] for atom in Fe_atoms]
        z = [atom.coordinates[2] for atom in Fe_atoms]
        
        position = [1 if (atom.label == "Fe(tet)" or atom.label == "Fe1") 
                    else 0 for atom in Fe_atoms]
        
        isAPB = [atom.isAPB for atom in Fe_atoms]
        unit_cell = [atom.ucn for atom in Fe_atoms]
        
        data = np.column_stack((x,y,z,position,isAPB,unit_cell))     
        np.savetxt(filename, data, fmt='%.5f')

    def rotate(self,alpha,beta,gamma):
        """ Rotate the crystal. """
        alpha = np.pi/180.0
        beta = np.pi/180.0
        gamma = np.pi/180.0
        
        rotation_x = np.array([[1.0,0.0,0.0],
                               [0.0, np.cos(alpha), -np.sin(alpha)],
                               [0.0, np.sin(alpha), np.cos(alpha)]])
        
        rotation_y = np.array([[np.cos(beta), 0.0, np.sin(beta)],
                               [0.0, 1.0, 0.0],
                               [-np.sin(beta), 0.0, np.cos(beta)]])
        
        rotation_z = np.array([[np.cos(gamma),-np.sin(gamma),0.0],
                               [np.sin(gamma), np.cos(gamma), 0.0],
                               [0.0, 0.0,1.0]])
        
        for atom in self.atoms:
            atom.coordinates -= self.radius
            atom.coordinates = rotation_x.dot(atom.coordinates)
            atom.coordinates = rotation_y.dot(atom.coordinates)
            atom.coordinates = rotation_z.dot(atom.coordinates)
            atom.coordinates += self.radius
  
