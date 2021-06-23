#include <iostream>
#include "randNumGenerator.hpp"
#include "Crystal.hpp"
#include "constants.h"
#include "Measurement.hpp"
#include "JSON.h"



#include <chrono>

struct retVals {
    double E_a, E_e, E_z, E_d;
};

retVals particle_configurations_test(std::string mode, std::string infile, std::string outfile, double c){
    std::string numf = "1";
    double center = c;
    std::string struct_file = infile;
//    if(APB){
//    struct_file ="/Users/tobiaskohler/PhD/APB-Paper/Monte-Carlo_Simulation/New_simulations_different_sizes/D11/D11_no_dip/5T_APB/D11_structure_APB_" + numf;
//    } else{
//        struct_file ="/Users/tobiaskohler/PhD/APB-Paper/Monte-Carlo_Simulation/New_simulations_different_sizes/D11/D11_no_dip/5T_noAPB/D11_structure_noAPB_" + numf;
//    }
    std::string spin_file = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/spins_" + numf + ".txt";
    Crystal crystal(struct_file, "brute_force", -21.0,-8.6,-28.1,-106.28,3.2e-25,0.0,0.0,-45.0, 0.1, center, 8.3965,8.3965,8.3965,0.03);
    
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
            if (crystal.atoms[i].y > center){
                if (crystal.atoms[i].position == 0){
                    crystal.atoms[i].spinx = -1.0;
                    crystal.atoms[i].spiny = 0.0;
                    crystal.atoms[i].spinz = 0.0;
                } else{
                    crystal.atoms[i].spinx = 1.0;
                    crystal.atoms[i].spiny = 0.0;
                    crystal.atoms[i].spinz = 0.0;
                }
            } else if (crystal.atoms[i].y <= center){
                if (crystal.atoms[i].position == 0){
                    crystal.atoms[i].spinx = 1.0;
                    crystal.atoms[i].spiny = 0.0;
                    crystal.atoms[i].spinz = 0.0;
                } else{
                    crystal.atoms[i].spinx = -1.0;
                    crystal.atoms[i].spiny = 0.0;
                    crystal.atoms[i].spinz = 0.0;
                }
            }
        }
    }
    
    if(mode == "aligned"){
        for(int i=0; i<crystal.atoms.size(); i++){
            if (crystal.atoms[i].position == 0){
                crystal.atoms[i].spinx = 1.0;
                crystal.atoms[i].spiny = 0.0;
                crystal.atoms[i].spinz = 0.0;
            } else{
                crystal.atoms[i].spinx = -1.0;
                crystal.atoms[i].spiny = 0.0;
                crystal.atoms[i].spinz = 0.0;
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
            crystal.atoms[i].spinx = u;
            crystal.atoms[i].spiny = v;
            crystal.atoms[i].spinz = w;
        }
        inFile.close();
        std::cout << "    Structure file reading successfull!" << std::endl;
    }
    
    for(int i =0; i<crystal.atoms.size(); i++){
        E_a += crystal.atoms[i].anisotropy();
        E_e += crystal.atoms[i].exchange();
        E_z += crystal.atoms[i].zeeman(5.0);
        E_d += crystal.atoms[i].dipole();
        sum_Sx += crystal.atoms[i].spinx;
    }
    E_total += E_a + E_e + E_z + E_d;
    std::cout << "E_a: " << E_a << std::endl;
    std::cout << "E_e: " << E_e << std::endl;
    std::cout << "E_z: " << E_z << std::endl;
    std::cout << "E_d: " << E_d << std::endl;
    std::cout << "E_total: " << E_total << std::endl;
    
    
    crystal.structure_snapshot(fname);
    return retVals{E_a, E_e, E_z, E_d};
}

void calc_dipole_field(){
    std::string numf = "1";
    double center = 6.0;
    std::string struct_file ="/Users/tobiaskohler/PhD/APB-Paper/Monte-Carlo_Simulation/New_simulations_different_sizes/D11/D11_no_dip/5T_noAPB/D11_structure_noAPB_" + numf;
    std::string spin_file = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/spins_" + numf + ".txt";
    Crystal crystal(struct_file, "brute_force", -21.0,-8.6,-28.1,-106.28,3.2e-25,0.0,0.0,-45.0, 0.1, center, 8.3965,8.3965,8.3965,0.03);
    
    std::string fname = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/90deg_domain_APB_D11.txt";
    std::ofstream dipole_field;
    dipole_field.open(fname, std::fstream::out);
    
    double H_dip[3];
    
    for(int i=0; i< crystal.atoms.size(); i++){
        crystal.atoms[i].dipole_field(H_dip);
        dipole_field << crystal.atoms[i].x << " " << crystal.atoms[i].y << " " << crystal.atoms[i].z << " " << H_dip[0] << " " << H_dip[1] << " " <<  H_dip[2] << std::endl;
    }
}

void relaxation_test(double size, double sigma, int steps){
    std::string path = "/Users/tobiaskohler/PhD/thesis/Simulations/tests/structures/";
    std::string struct_file = path + "D" + std::to_string((int)size) + "_structure_noAPB";
    std::string output_path = "/Users/tobiaskohler/PhD/thesis/Simulations/tests/relaxation_steps/";
    
    std::string sigma_str = std::to_string((int)(sigma)) + "d";
    sigma_str += std::to_string((int)(sigma*10.0-(int)(sigma)*10));
    sigma_str += std::to_string((int)(sigma*100.0-(int)(sigma*10.0)*10));
    
    
    double center = size/2.0+0.5;
    double temp = 0.01;
    Crystal crystal(struct_file, "None", -21.0,-8.6,-28.1,-106.28,3.2e-25,0.0,0.0,-45.0, 0.1, center, 8.3965,8.3965,8.3965,sigma);
    int numAtoms = crystal.atoms.size();
   
    for(int l=0; l<20; l++){
        crystal.reset_structure();
    std::ofstream relaxation;
    std::string output_file = output_path+ "relaxation_D" + std::to_string((int)size) + "_" + std::to_string(steps)+ "steps_" + sigma_str + "_" +std::to_string(l) + ".txt";
    relaxation.open(output_file, std::fstream::out);
    
    relaxation << "# MCS       Total trial moves      Sx comp.     Sy comp.    Sz comp." << std::endl;
    
    for(int j=0; j<steps; j++){
        for(int k=0; k<numAtoms; k++){
            crystal.atoms[rand0_crystalAtoms((int)crystal.atoms.size())].MonteCarloStep(5.0,temp);
        }
        double sum_sx = 0.0;
        double sum_sy = 0.0;
        double sum_sz = 0.0;
        for(int i =0; i<crystal.atoms.size(); i++){
            sum_sx += crystal.atoms[i].spinx;
            sum_sy += crystal.atoms[i].spiny;
            sum_sz += crystal.atoms[i].spinz;
        }
        relaxation << j << " " << j*numAtoms << " " << sum_sx <<" " << sum_sy << " " << sum_sz << std::endl;
    }
    }
}

void dipole_calcs(bool APB){
    std::ofstream dipole_calcs;
    std::string path = "/Users/tobiaskohler/PhD/APB-Paper/tests_structures_MC/dipole_energy_calculations/";
    if(APB){
        dipole_calcs.open(path+"dipole_calcs_APB.txt");
    } else{
        dipole_calcs.open(path+"dipole_calcs_noAPB.txt");
    }
    double E_d_aligned = 0.0;
    double E_d_domain = 0.0;
    
    dipole_calcs << "Diameter" << " " << "E_a_aligned" << " " << "E_a_domain" << " " << "E_e_aligned" << " " << "E_e_domain" << " " <<  "E_z_aligned" << " " << "E_z_domain" << " " << "E_d_aligned" << " " << "E_d_domain" <<  std::endl;
    
    for(int j=3; j<12; j++){
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
                                 path+out_aligned, c);
        retVals retVals_domain = particle_configurations_test("domain",
                                 path+f,
                                 path+out_domain, c);
//        dipole_calcs << j << " " << E_d_aligned << " " << E_d_domain << " " << E_a_aligned << " " << E_a_domain << " " <<  E_e_aligned << " " << E_e_domain << " " << E_z_aligned << " " << E_z_domain << std::endl;
        
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
    Crystal crystal(crystal_file, "None",
                    -21.0,-8.6,-28.1,-106.28,3.2e-25,0.0,0.0,-45.0, 0.1, 4.5,8.3965,8.3965,8.3965,sigma);
  
    double V[3];
    
    for(int i=0; i< 10000; i++){
        crystal.atoms[0].spinz = 1.0;
        crystal.atoms[0].spiny = 0.0;
        crystal.atoms[0].spinx = 0.0;
        
        crystal.atoms[0].angle(V);
        sigma_test << crystal.atoms[0].spinx << " " << crystal.atoms[0].spiny << " " << crystal.atoms[0].spinz << std::endl;
    }
}


int main(int argc, char* argv[]) {
    seed250(SEEED);
    
    if(argv[1] != NULL){
        std::string configfile = argv[1];
        std::ifstream ifs(configfile);
        nlohmann::json jf = nlohmann::json::parse(ifs);
        
    
        
        double macrocell_size;
        if(jf["dipole_interactions"]=="macrocell_method"){
            macrocell_size = jf["macrocell_size"];
        } else {
            macrocell_size = 0.0;
        }
    
        if(jf["Measurement"] == "M vs B"){
            run_MvsB(jf["output_dir"],
             jf["structure_file"],
             jf["dipole_interactions"],
             jf["steps"],
             jf["num_orientations"],
             jf["temperature"],
             jf["B_upper"],
             jf["B_lower"],
             jf["B_step"],
             jf["cooling_field"],
             jf["cooling_steps"],
             jf["start_temperature"],
             jf["temperature_step"],
             jf["FeTT"],
             jf["FeOO"],
             jf["FeTO"],
             jf["FeOO_APB"],
             jf["anisotropy constant"],
             macrocell_size,
             jf["lattice_a"],
             jf["lattice_b"],
             jf["lattice_c"],
             jf["particle center"],
             jf["sigma"]);
        } else if (jf["Measurement"] == "spin structure"){
            run_spinstructure(jf["dipole_interactions"],
                          jf["steps"],
                          jf["magnetic field"],
                          jf["temperature"],
                          jf["output_dir"],
                          jf["structure_file"],
                          jf["FeTT"],
                          jf["FeOO"],
                          jf["FeTO"],
                          jf["FeOO_APB"],
                          jf["anisotropy constant"],
                          jf["alpha"],
                          jf["beta"],
                          jf["gamma"],
                          macrocell_size,
                          jf["particle center"],
                          jf["lattice_a"],
                          jf["lattice_b"],
                          jf["lattice_c"],
                          jf["sigma"]);
        } else if(jf["Measurement"] == "M vs T"){
            run_MvsT(jf["output_dir"], 
                     jf["structure_file"],
                     jf["dipole_interactions"],
                     jf["steps"],
                     jf["averaging steps"],
                     jf["num_orientations"],
                     jf["measurement_field"],
                     jf["cooling_field"],
                     jf["TUpperLimit"],
                     jf["TLowerLimit"],
                     jf["TstepSize"],
                     jf["FeTT"],
                     jf["FeOO"],
                     jf["FeTO"],
                     jf["FeOO_APB"],
                     jf["anisotropy constant"],
                     jf["ZFC"],
                     jf["FC"],
                     macrocell_size,
                     jf["particle center"],
                     jf["lattice_a"],
                     jf["lattice_b"],
                     jf["lattice_c"],
                     jf["sigma"]);
        }
        else if (jf["Measurement"] == "test"){
            run_MvsB(jf["output_dir"],
             jf["structure_file"],
             "None",
             jf["steps"],
             jf["num_orientations"],
             jf["temperature"],
             jf["B_upper"],
             jf["B_lower"],
             jf["B_step"],
             jf["cooling_field"],
             jf["cooling_steps"],
             jf["start_temperature"],
             jf["temperature_step"],
             jf["FeTT"],
             jf["FeOO"],
             jf["FeTO"],
             jf["FeOO_APB"],
             jf["anisotropy constant"],
             macrocell_size,
              jf["lattice_a"],
               jf["lattice_b"],
              jf["lattice_c"],
           jf["particle center"],
             jf["sigma"]);
        }
    } else{
//        dipole_calcs(false);
//        dipole_calcs(true);
        
        
//        relaxation_test(6, 0.03, 5000);
//        relaxation_test(6, 0.1, 5000);
//        relaxation_test(6, 0.3, 5000);
//        
//        relaxation_test(8, 0.03, 5000);
//        relaxation_test(8, 0.1, 5000);
//        relaxation_test(8, 0.3, 5000);
        
//        relaxation_test(11, 0.03, 5000);
//        relaxation_test(11, 0.1, 5000);
//        relaxation_test(11, 0.3, 5000);
        
       // sigma_tests(0.4);
//        std::cout << "Please provide a configuration file in JSON format!" << std::endl;
//        std::cout << "For M vs B measurements the following arguments are needed:" << std::endl;
//        std::cout << "\t'structure_file'" << std::endl;
//        std::cout << "\t'dipole_interactions'" << std::endl;
//        std::cout << "\t'steps'" << std::endl;
//        std::cout << "\t'num_orientations'" << std::endl;
//        std::cout << "\t'temperature'" << std::endl;
//        std::cout << "\t'B_upper'" << std::endl;
//        std::cout << "\t'B_lower'" << std::endl;
//        std::cout << "\t'cooling_field'" << std::endl;
//        std::cout << "\t'cooling_steps'" << std::endl;
//        std::cout << "\t'start_temperature'" << std::endl;
//        std::cout << "\t'temperature_step'" << std::endl;
//        std::cout << "\t'FeTT'" << std::endl;
//        std::cout << "\t'FeOO'" << std::endl;
//        std::cout << "\t'FeTO'" << std::endl;
//        std::cout << "\t'FeOO_APB'" << std::endl;
//        std::cout << "\t'anisotropy constant'" << std::endl;
    }

    return 0;
}