//
//  Crystal2.hpp
//  MCS2
//
//  Created by Tobias Köhler on 06.01.21.
//  Copyright © 2021 Tobias Köhler. All rights reserved.
//

#ifndef Crystal_hpp
#define Crystal_hpp

#include <stdio.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include "Atom.hpp"
#include "randNumGenerator.hpp"
#include "constants.h"
#include "HelperStructs.hpp"
#include "Parsers.hpp"



class Crystal{
public:
    std::vector<Atom> atoms;
    LatticeParameters lattice_pars;
    double center;
    
    Crystal(std::string filename, DipoleInteractions dipoleInteractions,
            ExchangeConstants exchange_constants,
            double anisotropyConstant,
            LinalgVector angles,
            double macrocell_size, double center,
            LatticeParameters lattice_pars, double sigma);
    void rotateCrystal(LinalgVector, double);
    void randomOrientation();
    void alignAlongRandomVector();
    void structureSnapshot(std::string filename);
    void resetStructure();
    int outputStats();
    void setSigma(double sigma);
    void generateDipoleLists();
    void generateNeighbourLists();
    void initializeStructureFromFile(std::string filename,
                                  DipoleInteractions dipoleInteractions,
                                  ExchangeConstants exchange_constants,
                                  double anisotropyConstant);
    
    std::vector<Macrocell> macrocells;
    void generateMacrocells(double macrocell_size);
    void saveMacrocells(std::string filename);
    
};

#endif /* Crystal2_hpp */
