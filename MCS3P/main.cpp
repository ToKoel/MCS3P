#include <iostream>
#include "randNumGenerator.hpp"
#include "Crystal.hpp"
#include "constants.h"
#include "Measurement.hpp"
#include "JSON.h"
#include "testing.hpp"
#include <chrono>

int main(int argc, char* argv[]) {
    seed250(SEEED);
    
    if(argv[1] != NULL){
        // read configuration file
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
        dipole_calcs(false, 3, 11);
        dipole_calcs(true, 3, 11);
        
        
//        relaxation_test(6, 0.03, 8000, false);
//        relaxation_test(6, 0.1, 8000, false);
//        relaxation_test(6, 0.3, 8000, false);
//
//        relaxation_test(8, 0.03, 8000, false);
//        relaxation_test(8, 0.1, 8000, false);
//        relaxation_test(8, 0.3, 8000, false);
//
//        relaxation_test(11, 0.03, 8000, false);
//        relaxation_test(11, 0.1, 8000, false);
//        relaxation_test(11, 0.3, 8000, false);
//
//        relaxation_test(6, 0.03, 8000, true);
//        relaxation_test(6, 0.1, 8000, true);
//        relaxation_test(6, 0.3, 8000, true);
//
//        relaxation_test(8, 0.03, 8000, true);
//        relaxation_test(8, 0.1, 8000, true);
//        relaxation_test(8, 0.3, 8000, true);
//
//        relaxation_test(11, 0.03, 8000, true);
//        relaxation_test(11, 0.1, 8000, true);
//        relaxation_test(11, 0.3, 8000, true);
//
       // sigma_tests(0.4);
    }

    return 0;
}
