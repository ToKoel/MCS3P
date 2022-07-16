//
//  testing.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 27.06.21.
//

#include "testing.hpp"
MeasurementSettings measurementSettings;

retVals particle_configurations_test(std::string mode, std::string infile, std::string outfile, double c){
    std::string numf = "1";
    double center = c;
    std::string struct_file = infile;

    std::string spin_file = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/spins_" + numf + ".txt";
   // Crystal crystal(struct_file, "brute_force", -21.0,-8.6,-28.1,-106.28,3.25e-25,0.0,0.0,-45.0, 0.1, center, 8.3965,8.3965,8.3965,0.03);
    Crystal crystal(measurementSettings);
    
    // output file
    std::string fname = outfile;

    double E_a = 0.0;
    double E_e = 0.0;
    double E_z = 0.0;
    double E_d = 0.0;
    double E_total=0.0;
    double sum_Sx = 0.0;
    
    if(mode == "domain"){
        for(int i=0; i<crystal.atoms.size(); i++){
            if (crystal.atoms[i].positionVector.y > center){
                if (crystal.atoms[i].structuralPositionID == StructuralPositions::kOctahedral){
                    crystal.atoms[i].spinVector = {-1.0, 0.0, 0.0};
                } else{
                    crystal.atoms[i].spinVector = {1.0, 0.0, 0.0};
                }
            } else if (crystal.atoms[i].positionVector.y <= center){
                if (crystal.atoms[i].structuralPositionID == StructuralPositions::kOctahedral){
                    crystal.atoms[i].spinVector = {1.0, 0.0, 0.0};
                } else{
                    crystal.atoms[i].spinVector = {-1.0, 0.0, 0.0};
                }
            }
        }
    }
    
    if(mode == "aligned"){
        for(int i=0; i<crystal.atoms.size(); i++){
            if (crystal.atoms[i].structuralPositionID == StructuralPositions::kOctahedral){
                crystal.atoms[i].spinVector = {1.0, 0.0, 0.0};
            } else{
                crystal.atoms[i].spinVector = {-1.0, 0.0, 0.0};
            }
        }
    }
    
    if(mode == "calculated"){
        std::ifstream inFile;
        inFile.open(spin_file);
        std::string line;
        int num = 0;
        while (std::getline(inFile, line)){
                ++num;
        }
        inFile.close();
        std::cout << "    Number of lines in input file: " << num << std::endl;
        
        inFile.open(spin_file);
        for(int i =0; i<crystal.atoms.size(); i++){
            double u,v,w;
            inFile >> u >> v >> w ;
            crystal.atoms[i].spinVector = {u, v, w};
        }
        inFile.close();
        std::cout << "    Structure file reading successfull!" << std::endl;
    }
    
    for(int i =0; i<crystal.atoms.size(); i++){
        E_a += crystal.atoms[i].anisotropy();
        E_e += crystal.atoms[i].exchange();
        E_z += crystal.atoms[i].zeeman(5.0);
        E_d += crystal.atoms[i].dipole();
        sum_Sx += crystal.atoms[i].spinVector.x;
    }
    E_total += E_a + E_e + E_z + E_d;
    std::cout << "E_a: " << E_a << std::endl;
    std::cout << "E_e: " << E_e << std::endl;
    std::cout << "E_z: " << E_z << std::endl;
    std::cout << "E_d: " << E_d << std::endl;
    std::cout << "E_total: " << E_total << std::endl;
    
    
    crystal.structureSnapshot(fname);
    return retVals{E_a, E_e, E_z, E_d};
}

void calc_dipole_field(){
    std::string numf = "1";
    double center = 6.0;
    std::string struct_file ="/Users/tobiaskohler/PhD/APB-Paper/Monte-Carlo_Simulation/New_simulations_different_sizes/D11/D11_no_dip/5T_noAPB/D11_structure_noAPB_" + numf;
    std::string spin_file = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/spins_" + numf + ".txt";
    
    Crystal crystal(measurementSettings);
    
    std::string fname = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/90deg_domain_APB_D11.txt";
    std::ofstream dipole_field;
    dipole_field.open(fname, std::fstream::out);
    
    
    for(int i=0; i< crystal.atoms.size(); i++){
        LinalgVector H_dip = crystal.atoms[i].dipole_field();
        dipole_field << crystal.atoms[i].positionVector.x << " " << crystal.atoms[i].positionVector.y << " " << crystal.atoms[i].positionVector.z << " " << H_dip.x << " " << H_dip.y << " " <<  H_dip.z << std::endl;
    }
}

void relaxation_test(double size, double sigma, int steps, bool APB,
                     std::string output_path){
    std::string path = "/Users/tobiaskohler/PhD/thesis/Simulations/tests/structures/";
    std::string struct_file = path + "D" + std::to_string((int)size) + "_structure_";
    std::string A = " ";
    if(APB){
        struct_file += "APB";
        A = "_APB";
    } else{
        struct_file += "noAPB";
        A = "_noAPB";
    }
    
    std::cout << "Running relaxation test" << std::endl;
    
    std::string sigma_str = std::to_string((int)(sigma)) + "d";
    sigma_str += std::to_string((int)(sigma*10.0-(int)(sigma)*10));
    sigma_str += std::to_string((int)(sigma*100.0-(int)(sigma*10.0)*10));
    
    
    double center = size/2.0+0.5;
    double temp = 0.01;
  
    Crystal crystal(measurementSettings);
    
    int numAtoms = (int)crystal.atoms.size();
   
    for(int l=0; l<2; l++){
        crystal.resetStructure();
        std::ofstream relaxation;
        std::string output_file = output_path+ "relaxation_D" + std::to_string((int)size) + "_" + std::to_string(steps)+ "steps_" + sigma_str + "_" +std::to_string(l) + A + ".txt";
        relaxation.open(output_file, std::fstream::out);
    
        relaxation << "# MCS       Total trial moves      Sx comp.     Sy comp.    Sz comp." << std::endl;
    
        for(int j=0; j<steps; j++){
            for(int k=0; k<numAtoms; k++){
                crystal.atoms[rand0_crystalAtoms((int)crystal.atoms.size())].MonteCarloStep(5.0,temp);
            }
            LinalgVector sum{0.0,0.0,0.0};
            for(int i =0; i<crystal.atoms.size(); i++){
                sum += crystal.atoms[i].spinVector;
            }
            relaxation << j << " " << j*numAtoms << " " << sum.x <<" " << sum.y << " " << sum.z << std::endl;
        }
    }
}

void dipole_calcs(bool APB, int D_min, int D_max){
    std::ofstream dipole_calcs;
    std::string path = "/Users/tobiaskohler/PhD/thesis/Simulations/Energy_calculations/structures/";
    std::string path_out = "/Users/tobiaskohler/PhD/thesis/Simulations/Energy_calculations/calculations_wrong/";
    if(APB){
        dipole_calcs.open(path_out+"dipole_calcs_APB_1.txt");
    } else{
        dipole_calcs.open(path_out+"dipole_calcs_noAPB_1.txt");
    }
    
    dipole_calcs << "Diameter" << " " << "E_a_aligned" << " " << "E_a_domain" << " " << "E_e_aligned" << " " << "E_e_domain" << " " <<  "E_z_aligned" << " " << "E_z_domain" << " " << "E_d_aligned" << " " << "E_d_domain" <<  std::endl;
    
    for(int j=D_min; j<=D_max; j++){
    for(int i=1; i<21; i++){
        double c = j/2.0 + 0.5;
        std::string f = "D" + std::to_string(j);
        if(APB){
           f += "_structure_APB_" + std::to_string(i);
        } else{
           f += "_structure_noAPB_" + std::to_string(i);
        }
        std::string out_aligned = f + "_aligned.txt";
        std::string out_domain = f + "_domain.txt";
        retVals retVals_aligned = particle_configurations_test("aligned",
                                 path+f,
                                 path_out+out_aligned, c);
        retVals retVals_domain = particle_configurations_test("domain",
                                 path+f,
                                 path_out+out_domain, c);
        
        dipole_calcs << j << " " << retVals_aligned.E_a << " " << retVals_domain.E_a << " " << retVals_aligned.E_e << " " << retVals_domain.E_e << " " <<  retVals_aligned.E_z << " " << retVals_domain.E_z << " " << retVals_aligned.E_d << " " << retVals_domain.E_d << std::endl;
        }
    }
}

void sigma_tests(double sigma){
    std::ofstream sigma_test;
    std::string sigma_str = std::to_string((int)(sigma)) + "d";
    sigma_str += std::to_string((int)(sigma*10.0-(int)(sigma)*10));
    sigma_str += std::to_string((int)(sigma*100.0-(int)(sigma*10.0)*10));
    std::string outfile = "/Users/tobiaskohler/PhD/thesis/Simulations/ZFC_tests/sigma_tests/rand_vectors_gauss_" +sigma_str+".txt";

    sigma_test.open(outfile);
    std::string crystal_file = "/Users/tobiaskohler/PhD/thesis/Simulations/ZFC_tests/sigma_tests/D8_structure_noAPB_1";
    
    Crystal crystal(measurementSettings);
  
    
    for(int i=0; i< 10000; i++){
        crystal.atoms[0].spinVector = {1.0, 0.0, 0.0};
        
        crystal.atoms[0].angle();
        sigma_test << crystal.atoms[0].spinVector.x << " " << crystal.atoms[0].spinVector.y << " " << crystal.atoms[0].spinVector.z << std::endl;
    }
}
