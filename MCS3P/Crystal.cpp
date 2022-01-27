//
//  Crystal2.cpp
//  MCS2
//
//  Created by Tobias Köhler on 06.01.21.
//  Copyright © 2021 Tobias Köhler. All rights reserved.
//

#include "Crystal.hpp"


Crystal::Crystal(std::string filename, std::string dipole_interactions,
                 double FeTT, double FeOO, double FeTO, double FeOO_APB,
                 double anisotropyConstant, double alpha, double beta,
                 double gamma, double macrocell_size, double center,
                 double lattice_a, double lattice_b, double lattice_c, double sigma):
                 lattice_a(lattice_a), lattice_b(lattice_b), lattice_c(lattice_c), center(center){
    
    read_structure_from_file(filename, dipole_interactions,
                             FeTT,  FeOO,  FeTO,  FeOO_APB,
                             anisotropyConstant);
    set_sigma(sigma);
    rotateCrystal(alpha, beta, gamma, center);
    generate_neighbour_lists();

    if(dipole_interactions == "brute_force"){
        generate_dipole_lists();
    } else if(dipole_interactions == "macrocell_method"){
        generate_macrocells(macrocell_size);
    }
}

void Crystal::read_structure_from_file(std::string filename,
                                       std::string dipole_interactions,
                                       double FeTT, double FeOO, double FeTO, double FeOO_APB,
                                       double anisotropyConstant){

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<int> pos;
    std::vector<int> apb;
    
    std::ifstream inFile;
    inFile.open(filename);
    std::string line;
    int num = 0;
    while (std::getline(inFile, line)){
            ++num;
    }
    inFile.close();
    
    x.resize(num); y.resize(num); z.resize(num);
    pos.resize(num); apb.resize(num);
    
    inFile.open(filename);
    for(int i =0; i<num; i++){
        double xr,yr,zr;
        double posr,apbr,uc_xr,uc_yr,uc_zr;
        inFile >> xr >> yr >> zr >> posr >> apbr >> uc_xr >> uc_yr >> uc_zr ;
        x[i] = xr;
        y[i] = yr;
        z[i] = zr;
        pos[i] = posr;
        apb[i] = apbr;
    }
    inFile.close();
    
    int dip;
    if(dipole_interactions == "brute_force"){
        dip = 0;
    } else if(dipole_interactions == "macrocell_method"){
        dip = 1;
    } else {
        dip = 2;
    }
    
    for(int i=0; i<num; i++){
        atoms.push_back(Atom(dip,
                             x[i], y[i], z[i],
                             pos[i], apb[i],
                             FeTT,  FeOO,  FeTO, FeOO_APB,
                             anisotropyConstant));
    }
}

void Crystal::generate_neighbour_lists(){
    double xx, yy, zz;
    
    for(int i=0; i<atoms.size(); i++){
        atoms[i].neighboursAntiparallel.reserve(12);
        atoms[i].neighboursParallel.reserve(12);
        
        for(int j=0; j<atoms.size(); j++){
            xx = atoms[j].x - atoms[i].x;
            yy = atoms[j].y - atoms[i].y;
            zz = atoms[j].z - atoms[i].z;
    
            if((i!=j) and ((xx*xx + yy*yy + zz*zz) < 0.2505)){
                if(atoms[i].position == atoms[j].position){
                    atoms[i].neighboursParallel.push_back(&atoms[j]);
                }
                else{
                    atoms[i].neighboursAntiparallel.push_back(&atoms[j]);
                }
            }
        }
    }
}

void Crystal::generate_dipole_lists(){
    double x, y, z;
    double a_sq, b_sq, c_sq;
    double distance;
    double rx, ry, rz;
    
    a_sq = lattice_a*lattice_a;
    b_sq = lattice_b*lattice_b;
    c_sq = lattice_c*lattice_c;
    
    for(unsigned int i=0; i< atoms.size(); i++){
        x = atoms[i].x;
        y = atoms[i].y;
        z = atoms[i].z;
        for(unsigned int j=0; j< atoms.size(); j++){
            if(&atoms[i] != &atoms[j]){
                distance = std::sqrt((atoms[j].x - x)*(atoms[j].x - x)*a_sq+
                                     (atoms[j].y - y)*(atoms[j].y - y)*b_sq+
                                     (atoms[j].z - z)*(atoms[j].z - z)*c_sq);

                //unit vectors
                rx = (atoms[j].x - x)*lattice_a/distance;
                ry = (atoms[j].y - y)*lattice_b/distance;
                rz = (atoms[j].z - z)*lattice_c/distance;

                atoms[i].inv_distances_cubed.push_back(1.0/(std::pow((distance),3)));
                atoms[i].inv_distances_five.push_back(1.0/(std::pow((distance),5)));
                atoms[i].distVecX.push_back(rx);
                atoms[i].distVecY.push_back(ry);
                atoms[i].distVecZ.push_back(rz);

                atoms[i].magmag.push_back(MAGFE3*MAGFE3);

                atoms[i].allOtherAtomsInCrystal.push_back(&atoms[j]);
            } // end of distance calculation
        } // end of loop over other atoms
    } // end of loop over all atoms
}

int Crystal::outputStats(){
    int totalAtoms = (int)atoms.size();
    int apbAtoms = 0;
    int oct = 0;
    int tet = 0;

    for(int i =0; i< totalAtoms; i++){
        if(atoms[i].isApbAtom){
            apbAtoms++;
            if(atoms[i].position == 0){
                oct++;
            } else if(atoms[i].position == 1){
                tet++;
            }
        }
    }
    std::cout << "\nTotal number of Atoms: " << totalAtoms << std::endl;
    std::cout << "Number of APB atoms: " << apbAtoms << " (tet: "<< tet << "; oct: "<< oct << ")\n"<< std::endl;
    return totalAtoms;
}

void Crystal::structure_snapshot(std::string filename){
    std::ofstream structure;
    structure.open(filename, std::fstream::out);
    for(int i=0; i< atoms.size(); i++){
        structure << atoms[i].x << ", " << atoms[i].y << ", " << atoms[i].z <<", "<< atoms[i].spinx << ", " << atoms[i].spiny << ", " << atoms[i].spinz  << ", " << atoms[i].position << ", " << atoms[i].isApbAtom << "\n";
    }
    structure.close();
}

void Crystal::generate_macrocells(double macrocell_size){
    
    // find minimum position
    double x_min=atoms[0].x, y_min=atoms[0].y, z_min=atoms[0].z;
    for(int i=0; i<atoms.size(); i++){
        if(atoms[i].x < x_min){
            x_min = atoms[i].x;
        }
        if(atoms[i].y < y_min){
            y_min = atoms[i].y;
        }
        if(atoms[i].z < z_min){
            z_min = atoms[i].z;
        }
    }
    // find maximum position
    double x_max=atoms[0].x, y_max=atoms[0].y, z_max=atoms[0].z;
    for(int i=0; i<atoms.size(); i++){
        if(atoms[i].x > x_max){
            x_max = atoms[i].x;
        }
        if(atoms[i].y > y_max){
            y_max = atoms[i].y;
        }
        if(atoms[i].z > z_max){
            z_max = atoms[i].z;
        }
    }
    
    // determine number of macrocells needed
    int num_macrocells = 5;
    while(num_macrocells*macrocell_size <= x_max+0.2){
        num_macrocells++;
    }
    
    // initialize Macrocell instances
    for(int i=0; i<num_macrocells; i++){
        for(int j=0; j<num_macrocells; j++){
            for(int k=0; k<num_macrocells; k++){
                Macrocell macrocell(macrocell_size/2.0+(double)i*macrocell_size,
                                    macrocell_size/2.0+(double)j*macrocell_size,
                                    macrocell_size/2.0+(double)k*macrocell_size);
                macrocells.push_back(macrocell);
            }
        }
    }

    // assign atoms to macrocells
    double w = macrocell_size/2.0;
    for(int j=0; j<atoms.size(); j++){
        for(int i=0; i<macrocells.size(); i++){
            if(atoms[j].x < macrocells[i].center_x+w and atoms[j].x >= macrocells[i].center_x-w
               and atoms[j].y < macrocells[i].center_y+w and atoms[j].y >= macrocells[i].center_y-w
               and atoms[j].z < macrocells[i].center_z+w and atoms[j].z >= macrocells[i].center_z-w){
                macrocells[i].macrocell_atoms.push_back(&atoms[j]);
            }
        }
    }
    
    // flag and remove empty cells
    for(int i=0; i<macrocells.size(); i++){
        if(macrocells[i].macrocell_atoms.size()!=0){
            macrocells[i].isEmpty = false;
        } else {
            macrocells[i].isEmpty = true;
        }
    }
    macrocells.erase(std::remove_if(macrocells.begin(),
                                    macrocells.end(),
                                    [](Macrocell i){ return i.isEmpty==true;}),
                     macrocells.end());
    
    // Set pointers to macrocells for atoms
    for(int i=0; i<macrocells.size(); i++){
        for(int j=0; j<macrocells[i].macrocell_atoms.size(); j++){
            macrocells[i].macrocell_atoms[j]->macrocell_link = &macrocells[i];
        }
    }

    // shift macrocell center to center of mass and calculate effective volume
    for(int i=0; i<macrocells.size(); i++){
        double n = 0.0;
        macrocells[i].center_x = 0.0;
        macrocells[i].center_y = 0.0;
        macrocells[i].center_z = 0.0;
        // effective volume = N_atoms_per_macrocell*V_atom = N_atoms_per_macrocell * V_unit_cell/N_atoms_per_unit_cell
        macrocells[i].inv_effective_volume = 1.0/((macrocells[i].macrocell_atoms.size()*std::pow(8.3965,3)/24.0)*1e-30);
        for(int j=0; j<macrocells[i].macrocell_atoms.size(); j++){
            macrocells[i].center_x += macrocells[i].macrocell_atoms[j]->x;
            macrocells[i].center_y += macrocells[i].macrocell_atoms[j]->y;
            macrocells[i].center_z += macrocells[i].macrocell_atoms[j]->z;
            n+=1.0;
        }
        macrocells[i].center_x /= n;
        macrocells[i].center_y /= n;
        macrocells[i].center_z /= n;
    }
    
    // precalculate distances and distance vectors and initialize macrocell moment
    for(int i=0; i< macrocells.size(); i++){
        double x = macrocells[i].center_x;
        double y = macrocells[i].center_y;
        double z = macrocells[i].center_z;
        double distance = 0.0;
        double rx = 0.0, ry = 0.0, rz = 0.0;
        for(int j=0; j< macrocells.size(); j++){
            if(&macrocells[i] != &macrocells[j]){
                distance = std::sqrt((macrocells[j].center_x - x)*(macrocells[j].center_x - x)+
                                            (macrocells[j].center_y - y)*(macrocells[j].center_y - y)+
                                            (macrocells[j].center_z - z)*(macrocells[j].center_z - z));
              
                //unit vectors
                rx = (macrocells[j].center_x - x)/distance;
                ry = (macrocells[j].center_y - y)/distance;
                rz = (macrocells[j].center_z - z)/distance;

                macrocells[i].inv_distances_cubed.push_back(1.0/(std::pow((distance*8.3965),3)));
                macrocells[i].distVecX.push_back(rx);
                macrocells[i].distVecY.push_back(ry);
                macrocells[i].distVecZ.push_back(rz);

                macrocells[i].all_other_macrocells.push_back(&macrocells[j]);
            }
        }
        macrocells[i].total_moment[0] = 0.0;
        macrocells[i].total_moment[1] = 0.0;
        macrocells[i].total_moment[2] = 0.0;
        for(int k=0; k<macrocells[i].macrocell_atoms.size(); k++){
            macrocells[i].total_moment[0] += macrocells[i].macrocell_atoms[k]->spinx*MAGFE3;
            macrocells[i].total_moment[1] += macrocells[i].macrocell_atoms[k]->spiny*MAGFE3;
            macrocells[i].total_moment[2] += macrocells[i].macrocell_atoms[k]->spinz*MAGFE3;
        }
    }
}

void Crystal::save_macrocells(std::string filename){
    std::ofstream macrocells_centers;
    macrocells_centers.open(filename, std::fstream::out);
    for(int i=0; i< macrocells.size(); i++){
        macrocells_centers << macrocells[i].center_x << " "
                           << macrocells[i].center_y << " "
                           << macrocells[i].center_z << std::endl;
    }
    macrocells_centers.close();
}

void Crystal::reset_structure(){
    for(int i=0; i<atoms.size(); i++){
        double s[3];
        marsaglia(s);
        atoms[i].spinx = s[0];
        atoms[i].spiny = s[1];
        atoms[i].spinz = s[2];
    }
}

void Crystal::set_sigma(double sigma){
    for(int i=0; i< atoms.size(); i++){
        atoms[i].sigma = sigma;
    }
}

void Crystal::rotateCrystal(double alpha,double beta,double gamma, double center){
    
    double alphaRad = alpha * PI/180;
    double betaRad = beta * PI/180;
    double gammaRad = gamma * PI/180;
 
    double center_x = center;
    double center_y = center;
    double center_z = center;
    
    double x, y, z;
    double x2, y2, z2;
    double x3, y3, z3;


    for (unsigned int i = 0 ; i< atoms.size(); i++){
        // set particle into coordinate system origin
        x = atoms[i].x - center_x;
        y = atoms[i].y - center_y;
        z = atoms[i].z - center_z;

        // x-rotation with rotation matrix
        atoms[i].x = x*1.0 + y*0.0 + z*0.0;
        atoms[i].y = x*0.0 + y*std::cos(alphaRad) + z*-(std::sin(alphaRad));
        atoms[i].z = x*0.0 + y*std::sin(alphaRad) + z*std::cos(alphaRad);

        // temporary storing the new coordinates
        x2 = atoms[i].x;
        y2 = atoms[i].y;
        z2 = atoms[i].z;

        // y-rotation
        atoms[i].x = x2*std::cos(betaRad) + y2*0.0 + z2*std::sin(betaRad);
        atoms[i].y = x2*0.0 + y2*1.0 + z2*0.0;
        atoms[i].z = x2*-(std::sin(betaRad)) + y2*0.0 + z2*std::cos(betaRad);

        x3 = atoms[i].x;
        y3 = atoms[i].y;
        z3 = atoms[i].z;

        // z-rotation
        atoms[i].x = x3*std::cos(gammaRad) + y3*-(std::sin(gammaRad)) + z3*0.0;
        atoms[i].y = x3*std::sin(gammaRad) + y3*std::cos(gammaRad) + z3*0.0;
        atoms[i].z = x3*0.0 + y3*0.0 + z3*1.0;

        // particle gets set back to original position in space
        atoms[i].x = atoms[i].x + center_x;
        atoms[i].y = atoms[i].y + center_y;
        atoms[i].z = atoms[i].z + center_z;
    }
}

void Crystal::random_orientation(){
    double angles[3];
    rand0_360(angles);
    rotateCrystal(angles[0], angles[1], angles[2], center);
}

void Crystal::align_along_random_vector(){
    reset_structure();
    double random_magnetization_vector[3];
    marsaglia(random_magnetization_vector);
    for(int atom = 0; atom<(int)atoms.size(); atom++){
        if( atoms[atom].position == 0){
            atoms[atom].spinx = random_magnetization_vector[0];
            atoms[atom].spiny = random_magnetization_vector[1];
            atoms[atom].spinz = random_magnetization_vector[2];
        } else {
            atoms[atom].spinx = -random_magnetization_vector[0];
            atoms[atom].spiny = -random_magnetization_vector[1];
            atoms[atom].spinz = -random_magnetization_vector[2];
        }
    }
}

