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



class Crystal{
public:
    std::vector<Atom> atoms;
    // lattice parameters
    double lattice_a;
    double lattice_b;
    double lattice_c;
    
    Crystal(std::string filename, std::string dipole_interactions,
            double FeTT, double FeOO, double FeTO, double FeOO_APB,
            double anisotropyConstant, double alpha, double beta,
            double gamma, double macrocell_size, double center,
            double lattice_a, double lattice_b, double lattice_c, double sigma);
    void rotateCrystal(double, double, double, double);
    void structure_snapshot(std::string filename);
    void reset_structure();
    int outputStats();
    void set_sigma(double sigma);
    void generate_dipole_lists();
    void generate_neighbour_lists();
    void read_structure_from_file(std::string filename,
                                  std::string dipole_interactions,
                                  double FeTT, double FeOO, double FeTO, double FeOO_APB,
                                  double anisotropyConstant);
    
    std::vector<Macrocell> macrocells;
    void generate_macrocells(double macrocell_size);
    
};

#endif /* Crystal2_hpp */
