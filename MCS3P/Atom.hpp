#ifndef Atom_hpp
#define Atom_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "constants.h"
#include "randNumGenerator.hpp"
#include "omp.h"

class Atom;

class Macrocell{
public:
    double center_x;
    double center_y;
    double center_z;
    double total_moment[3];
    double previous_total_moment[3];
    double inv_effective_volume;
    bool isEmpty = true;
    std::vector<Atom*> macrocell_atoms;
    std::vector<Macrocell*> all_other_macrocells;
    
    std::vector<double> inv_distances_cubed;
    std::vector<double> distVecX;
    std::vector<double> distVecY;
    std::vector<double> distVecZ;
    
    Macrocell(double x,double y,double z);
    void update_total_moment();
    void get_demag_field(double*, double);
    void reset_total_moment();
};

class Atom {
    /*
     class to represent atomic magnetic moments with spatial position
     x,y and z, spin vector components spinx, spiny and spinz and
     other magnetic properties.
     */
public:
    double x = 0;   // atom coordinates
    double y = 0;
    double z = 0;
    double spinx = 0; // spin components
    double spiny = 0;
    double spinz = 0;
    double uc_x = 0.0; // unitCell_comp
    double uc_y = 0.0;
    double uc_z = 0.0;
    
    // exchange interaction constants
    double FeTT = 0.0;
    double FeOO = 0.0;
    double FeTO = 0.0;
    double FeOO_APB = 0.0;
    
    double anisotropyConstant = 0.0;
   
    int APB=0; // APB=1 --> APB position
    int position = 0;// Octahedral = 0 or tetrahedral = 1 position
    bool isApbAtom = false;
    bool isSurfaceAtom = false;
    int dipole_interactions; // 0 --> brute force, 1 --> macrocell method, 2 --> no dipole interactions
    bool macrocell_method;
    
    double sigma = 0.0; // defines opening angle of Gaussian cone
    double MagMagMu0 = MAGFE3*MAGFE3*MU0;
    
    Macrocell* macrocell_link = NULL;
    double H_demag[3] = {0.0,0.0,0.0};
    
    // Arrays containing pointers to Atom objects needed for the calculations
    std::vector<Atom*> neighboursAntiparallel;
    std::vector<Atom*> neighboursParallel;
    std::vector<Atom*> allOtherAtomsInCrystal;
    std::vector<double> inv_distances_cubed;
    std::vector<double> inv_distances_five;
    std::vector<double> distVecX;
    std::vector<double> distVecY;
    std::vector<double> distVecZ;
    std::vector<double> magmag;

    // initializer
    Atom(int dipole_interactions,
         double x, double y, double z,
         int position, int APB,
         double FeTT, double FeOO, double FeTO, double FeOO_APB,
         double anisotropyConstant);
    
    // energy functions
    double anisotropy();
    double exchange();
    double zeeman(double);
    double zeeman3D(double*);
    double dipole();
    
    void dipole_field(double*);
    
    // trial move functions
    void uniform(double *old_spin);
    void uniform_ziggurat(double *old_spin);
    void spin_flip(double*);
    void angle(double *old_spin);
    void hinzke_nowak(double *old_spin);
    void cattaneo_sun(double *old_spin);
    
    void MonteCarloStep(double Bx, double temperature);
};

#endif /* Atom_hpp */
