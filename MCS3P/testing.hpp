//
//  testing.hpp
//  MCS3P
//
//  Created by Tobias Köhler on 27.06.21.
//

#ifndef testing_hpp
#define testing_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Crystal.hpp"
#include "randNumGenerator.hpp"

struct retVals {
    double E_a, E_e, E_z, E_d;
};

void sigma_tests(double sigma);
void dipole_calcs(bool APB, int D_min, int D_max);
void relaxation_test(double size, double sigma, int steps, bool APB, std::string output_path);
void calc_dipole_field();
retVals particle_configurations_test(std::string mode, std::string infile, std::string outfile, double c);

#endif /* testing_hpp */
