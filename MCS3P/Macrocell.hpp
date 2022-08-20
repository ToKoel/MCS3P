////
////  Macrocell.hpp
////  MCS3P
////
////  Created by Tobias Köhler on 16.07.22.
////
//
//#ifndef Macrocell_hpp
//#define Macrocell_hpp
//
//#include <stdio.h>
//#include "Atom.hpp"
//
//class Macrocell{
//public:
//    double center_x;
//    double center_y;
//    double center_z;
//    double total_moment[3];
//    double previous_total_moment[3];
//    double inv_effective_volume;
//    bool isEmpty = true;
//    std::vector<Atom*> macrocell_atoms;
//    std::vector<Macrocell*> all_other_macrocells;
//    
//    std::vector<double> inv_distances_cubed;
//    std::vector<double> distVecX;
//    std::vector<double> distVecY;
//    std::vector<double> distVecZ;
//    
//    Macrocell(double x,double y,double z);
//    void update_total_moment();
//    void get_demag_field(LinalgVector&, double);
//    void reset_total_moment();
//};
//
//#endif /* Macrocell_hpp */
