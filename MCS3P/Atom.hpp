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
    LinalgVector tempSpinVector{0.0, 0.0, 0.0};
    LinalgVector unitCellVector{0.0, 0.0, 0.0};
    ExchangeConstants exchangeConstants{0.0, 0.0, 0.0, 0.0};
    double anisotropyConstant = 0.0;
    StructuralPositions structuralPositionID = StructuralPositions::kNone;
    DipoleInteractions dipoleInteractionHandling = DipoleInteractions::kNoInteractions;
    
    bool isApbAtom = false;
    bool isSurfaceAtom = false;
    
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

    Atom(MeasurementSettings,
         LinalgVector,
         StructuralPositions,
         bool);
    
    // energy functions
    double anisotropy();
    double exchange();
    double zeeman(double);
    double zeeman3D(LinalgVector B);
    double dipole();
    LinalgVector dipole_field();
    
    // trial move functions
    void uniform();
    void uniform_ziggurat();
    void spin_flip();
    void angle();
    void hinzke_nowak();
    void cattaneo_sun();
    
    void MonteCarloStep(double Bx, double temperature);
};

#endif /* Atom_hpp */
