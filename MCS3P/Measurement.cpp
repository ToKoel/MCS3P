#include "Measurement.hpp"


// function to generate field and temperature arrays
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    if(start < stop){
        for (T value = start; value < stop; value += step)
        values.push_back(value);
    } else{
        for (T value = start; value > stop; value -= step)
            values.push_back(value);
    }
    return values;
}

// helper function to convert doubles to string for filenames
std::string num_to_string(double number, int precision){
    std::string double_string = std::to_string((int)(number))+"d";
    double multiplier = 10.0;
    int A = 0;
    for(int i=0; i<precision; i++){
        A = (int)(number*multiplier - (int)(number*multiplier/10.0)*10);
        double_string += std::to_string(A);
        multiplier *= 10.0;
    }
    return double_string;
}

MvsTMeasurement::MvsTMeasurement(std::string dipoleInteractions,
                                 double steps,
                                 int averaging_steps,
                                 int numOrientations,
                                 double measurement_field,
                                 double cooling_field,
                                 double TUpperLimit,
                                 double TLowerLimit,
                                 double TstepSize,
                                 double macrocell_size):
dipoleInteractions(dipoleInteractions),
steps(steps),
averaging_steps(averaging_steps),
numOrientations(numOrientations),
measurement_field(measurement_field),
cooling_field(cooling_field),
TUpperLimit(TUpperLimit),
TLowerLimit(TLowerLimit),
TstepSize(TstepSize),
macrocell_size(macrocell_size){
}

void MvsTMeasurement::temperatureSweep(std::string output_dir, std::string structure_filename,
                                       double FeTT, double FeOO, double FeTO, double FeOO_APB,
                                       double anisotropyConstant,
                                       bool ZFC, bool FC,
                                       double center,
                                       double lattice_a, double lattice_b, double lattice_c,
                                       double sigma){
    
    std::string output_filename = output_dir;
    output_filename += structure_filename.substr(structure_filename.find_last_of("/")+1);
    
    output_filename += "_MvsT_sim_cF" + num_to_string(cooling_field, 2);
    output_filename += "T_mF" + num_to_string(measurement_field, 3);
    output_filename += "T_" + num_to_string(steps, 1) + "steps";
    output_filename += "_" + std::to_string((int)averaging_steps) + "avsteps";
    output_filename += "_" + std::to_string((int)numOrientations) + "or";
    output_filename += "_" + num_to_string(sigma, 2) + "sig";
    
    if(dipoleInteractions == "macrocell_method"){
        output_filename += "_0d" + std::to_string((int)(macrocell_size*100)) + "mcsize";
    }
    
    auto temperature_arr_down = arange<double>(TUpperLimit,TLowerLimit,TstepSize);
    auto temperature_arr_up = arange<double>(TLowerLimit,TUpperLimit,TstepSize);
    
    for(int i=0; i< temperature_arr_up.size();i++){
        mxZFC.push_back(0.0);
        myZFC.push_back(0.0);
        mzZFC.push_back(0.0);
    }
    for(int i=0; i< temperature_arr_down.size();i++){
        mxFC.push_back(0.0);
        myFC.push_back(0.0);
        mzFC.push_back(0.0);
    }
    
    ProgressBar bar;
    bar.set_bar_width(50);
    bar.fill_bar_progress_with("■");
    bar.fill_bar_remainder_with(" ");
    bar.update(0);
    
    Crystal crystal(structure_filename,
                    dipoleInteractions,
                    FeTT, FeOO, FeTO, FeOO_APB,
                    anisotropyConstant,
                    0.0, 0.0, 0.0, // orientation angles
                    macrocell_size, center,
                    lattice_a, lattice_b, lattice_c, sigma);

    for(int i = 0; i < numOrientations; i++){
        bar.update((double)(i)/(double)(numOrientations));
        bar.set_status_text("orientation " + std::to_string(i+1));
        
        // set random particle orientation
        crystal.random_orientation();
        // set spin structure fully aligned along random field
        crystal.align_along_random_vector();
        
        int totalNumAtoms = (int)crystal.atoms.size();
  
        if(ZFC){
            for(int j=0; j<temperature_arr_up.size(); j++){
            
                // relaxation steps
                for(int k = 0; k<(int)(steps*totalNumAtoms); k++){
                    crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field,
                                                                                    temperature_arr_up[j]);
                }
        
                double mx_temp = 0.0;
                double my_temp = 0.0;
                double mz_temp = 0.0;
            
                // averaging steps
                if(averaging_steps == 0){
                    for(int a=0; a<(int)crystal.atoms.size(); a++){
                        mx_temp += MAGFE3*crystal.atoms[a].spinx;
                        my_temp += MAGFE3*crystal.atoms[a].spiny;
                        mz_temp += MAGFE3*crystal.atoms[a].spinz;
                    }
                    mxZFC[j] += (mx_temp/(double)(totalNumAtoms))/(double)numOrientations;
                    myZFC[j] += (my_temp/(double)(totalNumAtoms))/(double)numOrientations;
                    mzZFC[j] += (mz_temp/(double)(totalNumAtoms))/(double)numOrientations;
                } else{
                    for(int m=0; m<averaging_steps*totalNumAtoms; m++){
                        crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field,
                                                                                        temperature_arr_up[j]);
                    
                        for(int a=0; a<(int)crystal.atoms.size(); a++){
                            mx_temp += MAGFE3*crystal.atoms[a].spinx;
                            my_temp += MAGFE3*crystal.atoms[a].spiny;
                            mz_temp += MAGFE3*crystal.atoms[a].spinz;
                        }
                    }
                    mxZFC[j] += (mx_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    myZFC[j] += (my_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    mzZFC[j] += (mz_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                } // end of averaging steps
            } // end of ZFC measurement
     
            // write parameters and results to file
            std::ofstream tSweepZfc;
            tSweepZfc.open (output_filename + "_ZFC.txt", std::fstream::out);
            tSweepZfc << "# structure_file: " << structure_filename << std::endl;
            tSweepZfc << "# dipole_interactions: " << dipoleInteractions << std::endl;
            tSweepZfc << "# steps: " << steps << std::endl;
            tSweepZfc << "# num_orientations: " << numOrientations  << std::endl;
            tSweepZfc << "# measurement_field: " << measurement_field << std::endl;
            tSweepZfc << "# cooling_field: " << cooling_field<< std::endl;
            tSweepZfc << "# TUpperLimit: " << TUpperLimit << std::endl;
            tSweepZfc << "# TLowerLimit: " << TLowerLimit << std::endl;
            tSweepZfc << "# TstepSize: " << TstepSize << std::endl;
            tSweepZfc << "# FeTT: " << FeTT << std::endl;
            tSweepZfc << "# FeTO: " << FeTO << std::endl;
            tSweepZfc << "# FeOO: " << FeOO << std::endl;
            tSweepZfc << "# FeOO_APB: " << FeOO_APB << std::endl;
            tSweepZfc << "# lattice parameters: " << lattice_a << ", " << lattice_b << ", " << lattice_c << std::endl;
            tSweepZfc << "# sigma: " << sigma << std::endl;
            for(int currentStep = 0; currentStep < temperature_arr_up.size(); currentStep++){
                tSweepZfc << temperature_arr_up[currentStep] << " " << mxZFC[currentStep] << " " << myZFC[currentStep] << " " << mzZFC[currentStep] << "\n";
            }
        } // end of ZFC
        
        if(FC){
            for(int j=0; j<temperature_arr_down.size(); j++){
          
                // relaxation steps
                for(int k = 0; k< (int)(steps*totalNumAtoms); k++){
                    crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field,
                                                                                    temperature_arr_down[j]);
                }
    
                double mx_temp_fc = 0.0;
                double my_temp_fc = 0.0;
                double mz_temp_fc = 0.0;

                // averaging steps
                if(averaging_steps == 0){
                    for(int a=0; a<(int)crystal.atoms.size(); a++){
                        mx_temp_fc += MAGFE3*crystal.atoms[a].spinx;
                        my_temp_fc += MAGFE3*crystal.atoms[a].spiny;
                        mz_temp_fc += MAGFE3*crystal.atoms[a].spinz;
                    }
                    mxFC[j] += (mx_temp_fc/(double)(totalNumAtoms))/(double)numOrientations;
                    myFC[j] += (my_temp_fc/(double)(totalNumAtoms))/(double)numOrientations;
                    mzFC[j] += (mz_temp_fc/(double)(totalNumAtoms))/(double)numOrientations;
                } else{
                    for(int m=0; m<averaging_steps*totalNumAtoms; m++){
                        crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field,
                                                                                        temperature_arr_down[j]);
                        for(int a=0; a<(int)crystal.atoms.size(); a++){
                            mx_temp_fc += MAGFE3*crystal.atoms[a].spinx;
                            my_temp_fc += MAGFE3*crystal.atoms[a].spiny;
                            mz_temp_fc += MAGFE3*crystal.atoms[a].spinz;
                        }
                    }
                    mxFC[j] += (mx_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    myFC[j] += (my_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    mzFC[j] += (mz_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                } // end of averaging steps
            } // end of FC measurement
           
            // write results and parameters to file
            std::ofstream tSweepFc;
            tSweepFc.open (output_filename + "_FC.txt", std::fstream::out);
            tSweepFc << "# structure_file: " << structure_filename << std::endl;
            tSweepFc << "# dipole_interactions: " << dipoleInteractions << std::endl;
            tSweepFc << "# steps: " << steps << std::endl;
            tSweepFc << "# num_orientations: " << numOrientations  << std::endl;
            tSweepFc << "# measurement_field: " << measurement_field << std::endl;
            tSweepFc << "# cooling_field: " << cooling_field<< std::endl;
            tSweepFc << "# TUpperLimit: " << TUpperLimit << std::endl;
            tSweepFc << "# TLowerLimit: " << TLowerLimit << std::endl;
            tSweepFc << "# TstepSize: " << TstepSize << std::endl;
            tSweepFc << "# FeTT: " << FeTT << std::endl;
            tSweepFc << "# FeTO: " << FeTO << std::endl;
            tSweepFc << "# FeOO: " << FeOO << std::endl;
            tSweepFc << "# FeOO_APB: " << FeOO_APB << std::endl;
            tSweepFc << "# sigma: " << sigma << std::endl;
            tSweepFc << "# lattice parameters: " << lattice_a << ", " << lattice_b << ", " << lattice_c << std::endl;
            tSweepFc << "# " << std::endl;
            tSweepFc << "# field (T)      M_x(B)      M_y(B)      M_z(B)" << std::endl;
            for(int currentStep = 0; currentStep < temperature_arr_down.size(); currentStep++){
                tSweepFc << temperature_arr_down[currentStep] << " " << mxFC[currentStep] << " " << myFC[currentStep] << " " << mzFC[currentStep] << "\n";
            }
        } // end of FC
    } // end of orientation
} // end of temperature sweep

void run_MvsT(std::string output_dir, std::string structure_filename,std::string dipoleInteractions,
              double steps, int averaging_steps, int numOrientations, double measurement_field,
              double cooling_field, double TUpperLimit, double TLowerLimit, double TstepSize,
              double FeTT, double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant,
              bool ZFC, bool FC, double macrocell_size, double center,double lattice_a,double lattice_b,
              double lattice_c, double sigma){
    
    MvsTMeasurement MvsTMeasurement(dipoleInteractions,
                                    steps,
                                    averaging_steps,
                                    numOrientations,
                                    measurement_field,
                                    cooling_field,
                                    TUpperLimit,
                                    TLowerLimit,
                                    TstepSize, macrocell_size);
    MvsTMeasurement.temperatureSweep(output_dir,structure_filename,
                                      FeTT,  FeOO,  FeTO,  FeOO_APB,
                                      anisotropyConstant, ZFC, FC, center,
                                     lattice_a, lattice_b, lattice_c,
                                     sigma);
}


MvsBMeasurement::MvsBMeasurement(std::string dipoleInteractions,
                                 int steps,
                                 int numOrientations,
                                 double temperature,
                                 double BUpperLimit,
                                 double BLowerLimit,
                                 double BstepSize,
                                 double startTemp,
                                 double tempStep,
                                 double coolingField,
                                 int coolingSteps, double macrocell_size):
dipoleInteractions(dipoleInteractions),
steps(steps),
numOrientations(numOrientations),
temperature(temperature),
BUpperLimit(BUpperLimit), BLowerLimit(BLowerLimit), BstepSize(BstepSize),
startTemp(startTemp), tempStep(tempStep),
coolingField(coolingField), coolingSteps(coolingSteps),
macrocell_size(macrocell_size) {
    
    // generate field arrays
    BnumberOfSteps = (BUpperLimit/BstepSize)+2*((BUpperLimit-BLowerLimit)/BstepSize);
    int BstepsUpStart = BUpperLimit/BstepSize;
    int BstepsDown = (BUpperLimit - BLowerLimit)/BstepSize;
    for(int i = 0; i<BstepsUpStart; i++){
        magneticField.push_back(0.0 + i*BstepSize);
    }
    for(int i = 0; i < BstepsDown; i++){
        magneticField.push_back(BUpperLimit - i*BstepSize);
    }
    for(int i = 0; i < BstepsDown; i++){
        magneticField.push_back(BLowerLimit + i*BstepSize);
    }
    magneticField.push_back(BUpperLimit);
    
    // generate magnetization vectors for each field step
    for(int i = 0; i <= BnumberOfSteps; i++){
        mX.push_back(0.0);
        mY.push_back(0.0);
        mZ.push_back(0.0);
    }
}

void MvsBMeasurement::fieldSweep(std::string output_dir,
                                 std::string structure_filename,
                                 double FeTT,double FeOO, double FeTO, double FeOO_APB,
                                 double anisotropyConstant,
                                 double lattice_a, double lattice_b, double lattice_c,
                                 double center, double sigma){
    
    std::string output_filename = structure_filename.substr(structure_filename.find_last_of("/")+1);
    output_filename += "_Hysteresis_sim_"+num_to_string(temperature, 1);
    output_filename += "K_"+std::to_string((int)(steps))+"steps";
    output_filename += "_"+std::to_string((int)(numOrientations))+"or";
    output_filename += "_"+num_to_string(coolingField, 2)+"T";
    output_filename += "_"+num_to_string(coolingSteps, 2)+"cs";
    output_filename += "_sT"+num_to_string(startTemp, 1)+"K";
    output_filename += "_dip"+dipoleInteractions;
    if(dipoleInteractions == "macrocell_method"){
        output_filename += "_0d" + std::to_string((int)(macrocell_size*100)) + "mcsize";
    }
    
    std::cout << "\nOutput filename: " << output_filename << std::endl;
    std::cout << "\n ------- starting M vs B measurement -------\n" << std::endl;
    
    output_filename = output_dir + output_filename;
    
    for(int i=0; i <= BnumberOfSteps; i++){
        mX[i] = 0.0;
        mY[i] = 0.0;
        mZ[i] = 0.0;
    }
    
    ProgressBar bar;
    bar.set_bar_width(50);
    bar.fill_bar_progress_with("■");
    bar.fill_bar_remainder_with(" ");
    bar.update(0);
    
    for(int i = 0; i < numOrientations; i++){
        bar.update((double)(i)/(double)(numOrientations));
        bar.set_status_text("orientation " + std::to_string(i+1));
        
        double angles[3];
        rand0_360(angles);
    
        Crystal crystal(structure_filename, dipoleInteractions,
                        FeTT, FeOO, FeTO, FeOO_APB,
                        anisotropyConstant,
                        angles[0], angles[1], angles[2],
                        macrocell_size, center,
                        lattice_a, lattice_b, lattice_c,
                        sigma);
        
        int totalNumAtoms = int(crystal.atoms.size());
   
        // (field) cooling
        auto temperature_arr = arange<double>(startTemp,temperature,tempStep);
        for(int i=0; i<temperature_arr.size();i++){
            for(int j=0; j<coolingSteps * totalNumAtoms; j++){
                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(coolingField,
                                                                                temperature_arr[i]);
            }
        }
 
        // field sweep
        for(int j=0; j<=BnumberOfSteps; j++){
            for(int k = 0; k<steps*totalNumAtoms; k++){
                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(magneticField[j],
                                                                                temperature);
            }
            for(unsigned int l=0; l<=crystal.atoms.size(); l++){
                mX[j] += MAGFE3*crystal.atoms[l].spinx;
                mY[j] += MAGFE3*crystal.atoms[l].spiny;
                mZ[j] += MAGFE3*crystal.atoms[l].spinz;
            }
        }
        
        std::ofstream bSweep;
        bSweep.open (output_filename + ".txt", std::fstream::out);
        bSweep << "# structure_file: " << structure_filename << std::endl;
        bSweep << "# dipole_interactions: " << dipoleInteractions << std::endl;
        bSweep << "# steps: " << steps << std::endl;
        bSweep << "# num_orientations: " << numOrientations  << std::endl;
        bSweep << "# temperature: " << temperature << std::endl;
        bSweep << "# B_upper: " << BUpperLimit<< std::endl;
        bSweep << "# B_lower: " << BLowerLimit << std::endl;
        bSweep << "# B_step: " << BstepSize << std::endl;
        bSweep << "# cooling_field: " << coolingField << std::endl;
        bSweep << "# start_temperature: " << startTemp << std::endl;
        bSweep << "# temperature_step: " << tempStep << std::endl;
        bSweep << "# FeTT: " << FeTT << std::endl;
        bSweep << "# FeTO: " << FeTO << std::endl;
        bSweep << "# FeOO: " << FeOO << std::endl;
        bSweep << "# FeOO_APB: " << FeOO_APB << std::endl;
        bSweep << "# lattice parameters: " << lattice_a << ", " << lattice_b << ", " << lattice_c << std::endl;
        bSweep << "# total number of iron atoms: " << totalNumAtoms << std::endl;
        bSweep << "# " << std::endl;
        bSweep << "# field (T)      M_x(B)      M_y(B)      M_z(B)" << std::endl;
        for(int currentStep = 0; currentStep <= BnumberOfSteps; currentStep++){
            bSweep << magneticField[currentStep] << " " << mX[currentStep] << " " << mY[currentStep] << " " << mZ[currentStep] << "\n";
        }
    } // end of field sweep for one particle orientation
} // end of field sweep

void run_MvsB(std::string output_dir,
              std::string structure_filename,
              std::string dipoleInteractions,
              int steps, int numOrientations,
              double temperature,
              double BUpperLimit, double BLowerLimit, double BstepSize,
              double coolingField, int coolingSteps, double startTemp, double tempStep,
              double FeTT, double FeOO, double FeTO, double FeOO_APB,
              double anisotropyConstant,
              double macrocell_size,
              double lattice_a, double lattice_b, double lattice_c,
              double center, double sigma){
    MvsBMeasurement MvsBMeasurement(dipoleInteractions,
                                    steps,
                                    numOrientations,
                                    temperature,
                                    BUpperLimit,
                                    BLowerLimit,
                                    BstepSize,
                                    startTemp,
                                    tempStep,
                                    coolingField,
                                    coolingSteps,
                                    macrocell_size);
    MvsBMeasurement.fieldSweep(output_dir, structure_filename,
                               FeTT,  FeOO,  FeTO,  FeOO_APB,
                               anisotropyConstant,
                               lattice_a,  lattice_b,  lattice_c,
                               center, sigma);
}

spinStructure::spinStructure(std::string dipoleInteractions,
                             int steps,
                             double magneticField,
                             double temperature,
                             double macrocell_size,
                             std::string output_dir,
                             std::string structure_filename):
                             dipoleInteractions(dipoleInteractions),
                             steps(steps),
                             magneticField(magneticField),
                             temperature(temperature),
                             macrocell_size(macrocell_size),
                             output_dir(output_dir),
                             structure_filename(structure_filename){
};

void spinStructure::spinStructureMeasurement(double FeTT,double FeOO, double FeTO, double FeOO_APB,
                                             double anisotropyConstant,
                                             double alpha, double beta, double gamma,
                                             double center,
                                             double lattice_a, double lattice_b, double lattice_c,
                                             double sigma){
    Crystal crystal(structure_filename, dipoleInteractions,
                    FeTT, FeOO, FeTO, FeOO_APB,
                    anisotropyConstant,
                    alpha, beta, gamma,
                    macrocell_size, center,
                    lattice_a, lattice_b, lattice_c,
                    sigma);

    int totalNumAtoms = int(crystal.atoms.size());
    
    ProgressBar bar2;
    bar2.set_bar_width(50);
    bar2.fill_bar_progress_with("■");
    bar2.fill_bar_remainder_with(" ");
    bar2.update(0);
    
    //simulated annealing
    /*
    std::cout<< "cooling" << std::endl;
    bar2.set_status_text("starting cooling");
    double start_T = 300.0;
    while(start_T > temperature){
        for(int i=0; i<steps*totalNumAtoms; i++){
            crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(magneticField, start_T);
        }
        start_T -= 10.0;
        std::cout << start_T << std::endl;
        bar2.set_status_text("temperature: " + std::to_string((int)start_T));
    }
    bar2.set_status_text("cooling finished");
    */
    
    ProgressBar bar;
    bar.set_bar_width(50);
    bar.fill_bar_progress_with("■");
    bar.fill_bar_remainder_with(" ");
    bar.update(0);
    
    // measurement
    {
        for(int i =0; i < steps*totalNumAtoms; i++){
            crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(magneticField, temperature);
            if(i % 100000 == 0){
                bar.update((double)(i)/(double)(steps*totalNumAtoms));
                bar.set_status_text("trial move " + std::to_string(i));
            }
        }

            
        std::string output_filename = output_dir;
        output_filename += structure_filename.substr(structure_filename.find_last_of("/")+1);
        output_filename += "_spin_structure_"+std::to_string((int)(temperature));
        output_filename += "K_"+std::to_string((int)(steps))+"MCS";
        output_filename += "_"+ std::to_string((int)magneticField) +"T";
        output_filename += "_dip"+dipoleInteractions;
        if(dipoleInteractions == "macrocell_method"){
            output_filename += "_0d" + std::to_string((int)(macrocell_size*100)) + "mcsize";
        }
        
        std::ofstream structure;
        structure.open(output_filename + ".txt", std::fstream::out);
        
        structure << "# structure_file: " << structure_filename << std::endl;
        structure << "# dipole_interactions: " << dipoleInteractions << std::endl;
        if(dipoleInteractions == "macrocell_method"){
            structure << "# macrocell_size: " << macrocell_size << std::endl;
        }
        structure << "# steps: " << steps << std::endl;
        structure << "# Number of atoms: " << totalNumAtoms << std::endl;
        structure << "# temperature: " << temperature << std::endl;
        structure << "# field: " << magneticField << std::endl;
        structure << "# FeTT: " << FeTT << std::endl;
        structure << "# FeTO: " << FeTO << std::endl;
        structure << "# FeOO: " << FeOO << std::endl;
        structure << "# FeOO_APB: " << FeOO_APB << std::endl;
        structure << "# anisotropy constant: " << anisotropyConstant << std::endl;
        structure << "# lattice parameters: " << lattice_a << ", " << lattice_b << ", " << lattice_c << std::endl;
        structure << "# Orientation: " << alpha << ", " << beta << ", " << gamma << std::endl;
        structure << "# " << std::endl;
        structure << "# x      y        z      spinx    spiny    spinz  pos  APB" << std::endl;
        
        for(int i=0; i< crystal.atoms.size(); i++){
            structure << crystal.atoms[i].x << ", "
                      << crystal.atoms[i].y << ", "
                      << crystal.atoms[i].z << ", "
                      << crystal.atoms[i].spinx << ", "
                      << crystal.atoms[i].spiny << ", "
                      << crystal.atoms[i].spinz << ", "
                      << crystal.atoms[i].position << ", "
                      << crystal.atoms[i].isApbAtom << "\n";
        }
    } // end of spin structure simulation
}

void run_spinstructure(std::string dipoleInteractions,
                       int steps, double magneticField,
                       double temperature,
                       std::string output_dir,
                       std::string structure_filename,
                       double FeTT,double FeOO, double FeTO, double FeOO_APB,
                       double anisotropyConstant,
                       double alpha,double beta, double gamma,
                       double macrocell_size, double center,
                       double lattice_a, double lattice_b, double lattice_c,
                       double sigma){
    spinStructure spinStructure(dipoleInteractions,
                                steps,
                                magneticField,
                                temperature,
                                macrocell_size,
                                output_dir,
                                structure_filename
                                );
    spinStructure.spinStructureMeasurement(FeTT, FeOO, FeTO, FeOO_APB,
                                           anisotropyConstant, alpha,
                                           beta, gamma, center,  lattice_a,
                                           lattice_b, lattice_c, sigma);
}
