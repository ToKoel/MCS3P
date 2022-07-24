#ifndef Atom_hpp
#define Atom_hpp

#include <vector>
#include <iostream>
#include <cmath>
#include "constants.hpp"
#include "HelperStructs.hpp"
#include "randNumGenerator.hpp"
#include "omp.h"

class Macrocell;

class Atom {
public:
    utility::LinalgVector positionVector{0.0, 0.0, 0.0};
    utility::LinalgVector spinVector{0.0, 0.0, 0.0};
    utility::LinalgVector tempSpinVector{0.0, 0.0, 0.0};
    utility::LinalgVector unitCellVector{0.0, 0.0, 0.0};
    utility::ExchangeConstants exchangeConstants{0.0, 0.0, 0.0, 0.0};
    double anisotropyConstant = 0.0;
    utility::StructuralPositions structuralPositionID = utility::StructuralPositions::kNone;
    utility::DipoleInteractions dipoleInteractionHandling = utility::DipoleInteractions::kNoInteractions;
    
    bool isApbAtom = false;
    bool isSurfaceAtom = false;
    
    double sigma = 0.0; // defines opening angle of Gaussian cone
    double MagMagMu0 = MAGFE3*MAGFE3*MU0;
    
    Macrocell* macrocell_link = nullptr;
    utility::LinalgVector H_demag = {0.0,0.0,0.0};
    
    // Arrays containing pointers to Atom objects needed for the calculations
    std::vector<Atom*> neighboursAntiparallel;
    std::vector<Atom*> neighboursParallel;
    std::vector<Atom*> allOtherAtomsInCrystal;
    std::vector<double> inv_distances_cubed;
    std::vector<double> inv_distances_five;
    std::vector<utility::LinalgVector> distanceVectors;
    std::vector<double> magmag;

    Atom(const utility::MeasurementSettings& measurementSettings, utility::LinalgVector positionVector,
         utility::StructuralPositions structuralPositionID, bool isAPB);
    
    // energy functions
    double const anisotropy();
    double const exchange();
    double const zeeman(const double magneticFieldXcomponent);
    double const zeeman3D(const utility::LinalgVector magneticFieldVector);
    double const dipole();
    utility::LinalgVector const dipole_field();
    
    // trial move functions
    void uniform();
    void uniform_ziggurat();
    void spin_flip();
    void angle();
    void hinzke_nowak();
    void cattaneo_sun();
    
    void MonteCarloStep(utility::Environment environment);
};

#endif /* Atom_hpp */
