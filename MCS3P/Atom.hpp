#ifndef Atom_hpp
#define Atom_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "constants.h"
#include "HelperStructs.hpp"
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
    void get_demag_field(LinalgVector&, double);
    void reset_total_moment();
};

class Atom {
    /*
     class to represent atomic magnetic moments with spatial position
     x,y and z, spin vector components spinx, spiny and spinz and
     other magnetic properties.
     */
public:
    LinalgVector positionVector{0.0, 0.0, 0.0};
    LinalgVector spinVector{0.0, 0.0, 0.0};
    LinalgVector unitCellVector{0.0, 0.0, 0.0};
    ExchangeConstants exchangeConstants{0.0, 0.0, 0.0, 0.0};
    double anisotropyConstant = 0.0;
   
    StructuralPositions structuralPositionID = StructuralPositions::kOctahedral;
    DipoleInteractions dipoleInteractionHandling = DipoleInteractions::kNoInteractions;
    
    bool isApbAtom = false;
    bool isSurfaceAtom = false;
    bool macrocell_method = false;
    
    double sigma = 0.0; // defines opening angle of Gaussian cone
    double MagMagMu0 = MAGFE3*MAGFE3*MU0;
    
    Macrocell* macrocell_link = nullptr;
    LinalgVector H_demag = {0.0,0.0,0.0};
    
    // Arrays containing pointers to Atom objects needed for the calculations
    std::vector<Atom*> neighboursAntiparallel;
    std::vector<Atom*> neighboursParallel;
    std::vector<Atom*> allOtherAtomsInCrystal;
    std::vector<double> inv_distances_cubed;
    std::vector<double> inv_distances_five;
    std::vector<LinalgVector> distanceVectors;
    std::vector<double> magmag;

    Atom(DipoleInteractions dipoleInteractionHandling,
         LinalgVector positionVector,
         StructuralPositions structuralPositionID,
         bool isAPB,
         ExchangeConstants exchangeConstants,
         double anisotropyConstant);
    
    // energy functions
    double anisotropy();
    double exchange();
    double zeeman(double);
    double zeeman3D(LinalgVector B);
    double dipole();
    void dipole_field(double*);
    
    // trial move functions
    void uniform(LinalgVector& oldSpin);
    void uniform_ziggurat(LinalgVector& oldSpin);
    void spin_flip(LinalgVector& oldSpin);
    void angle(LinalgVector& oldSpin);
    void hinzke_nowak(LinalgVector& oldSpin);
    void cattaneo_sun(LinalgVector& oldSpin);
    
    void MonteCarloStep(double Bx, double temperature);
};

#endif /* Atom_hpp */
