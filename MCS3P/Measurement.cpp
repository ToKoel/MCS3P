#include "Measurement.hpp"

class Timer{
public:
    Timer(){
        m_StartTimepoint = std::chrono::high_resolution_clock::now();
    }
    
    ~Timer(){
        Stop();
    }
    
    void Stop(){
        auto endTimepoint = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
        
        auto duration = end-start;
        
        std::cout << duration/1000000.0 << "s\n";
    }
private:
    std::chrono::time_point< std::chrono::high_resolution_clock> m_StartTimepoint;
};


Measurement::Measurement(int particleSize, bool antiphaseBoundary, double surfaceLayerThickness, bool dipoleInteractions,bool vacancies, bool cutOff, int steps): particleSize(particleSize), antiphaseBoundary(antiphaseBoundary), surfaceLayerThickness(surfaceLayerThickness), dipoleInteractions(dipoleInteractions),vacancies(vacancies), cutOff(cutOff), steps(steps){
}


MvsTMeasurement::MvsTMeasurement(std::string dipoleInteractions,double steps,int averaging_steps, int numOrientations, double measurement_field, double cooling_field, double TUpperLimit, double TLowerLimit, double TstepSize, double macrocell_size):dipoleInteractions(dipoleInteractions),steps(steps), averaging_steps(averaging_steps), numOrientations(numOrientations), measurement_field(measurement_field),cooling_field(cooling_field), TUpperLimit(TUpperLimit), TLowerLimit(TLowerLimit), TstepSize(TstepSize), macrocell_size(macrocell_size){
    
    TnumberOfSteps = (TUpperLimit-TLowerLimit)/TstepSize;
    for (int i = 0; i < TnumberOfSteps; i++){
        Tdown.push_back(TUpperLimit-(TstepSize*i));
        Tup.push_back(TLowerLimit+(TstepSize*i));
    }
}

void MvsTMeasurement::temperatureSweep(std::string output_dir, std::string structure_filename,
                                       double FeTT, double FeOO, double FeTO, double FeOO_APB,
                                       double anisotropyConstant,bool ZFC, bool FC, double center,double lattice_a, double lattice_b, double lattice_c, double sigma){
    std::string output_filename = output_dir + structure_filename.substr(structure_filename.find_last_of("/")+1);
    
    std::string cF = std::to_string((int)(cooling_field))+"d";
    cF += std::to_string((int)(cooling_field*10.0-(int)(cooling_field)*10));
    cF += std::to_string((int)(cooling_field*100.0-(int)(cooling_field*10.0)*10));
    cF += std::to_string((int)(cooling_field*1000.0-(int)(cooling_field*100.0)*10));
    
    std::string mF = std::to_string((int)(measurement_field))+"d";
    mF += std::to_string((int)(measurement_field*10.0-(int)(measurement_field)*10));
    mF += std::to_string((int)(measurement_field*100.0-(int)(measurement_field*10.0)*10));
    mF += std::to_string((int)(measurement_field*1000.0-(int)(measurement_field*100.0)*10));
    
    std::string steps_str = std::to_string((int)(steps)) + "d";
    steps_str += std::to_string((int)(steps*10.0-(int)(steps)*10));
    steps_str += std::to_string((int)(steps*100.0-(int)(steps*10.0)*10));
    
    std::string sigma_str = std::to_string((int)(sigma)) + "d";
    sigma_str += std::to_string((int)(sigma*10.0-(int)(sigma)*10));
    sigma_str += std::to_string((int)(sigma*100.0-(int)(sigma*10.0)*10));
 
    output_filename += "_MvsT_sim_cF"+cF+"T_mF"+mF;
    output_filename += "T_"+steps_str+"steps";
    output_filename += "_"+std::to_string((int)averaging_steps)+"avsteps";
    output_filename += "_"+std::to_string((int)(numOrientations))+"or";
    output_filename += "_"+sigma_str+"sig";
    
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
    
    Crystal crystal(structure_filename, dipoleInteractions,
                    FeTT, FeOO, FeTO, FeOO_APB, anisotropyConstant,
                    0.0, 0.0, 0.0,macrocell_size, center,lattice_a, lattice_b, lattice_c, sigma);
    
    //std::ofstream orientations_log;
    //orientations_log.open(output_filename + "_orientations.txt", std::fstream::out);
    
//    for(int a=0; a<(int)(3000*(int)crystal.atoms.size()); a++){
//        crystal.atoms[rand0_crystalAtoms((int)crystal.atoms.size())].MonteCarloStep(cooling_field, temperature_arr_up[0]);
//    }

    for(int i = 0; i < numOrientations; i++){
        bar.update((double)(i)/(double)(numOrientations));
        bar.set_status_text("orientation " + std::to_string(i+1));
        
        double angles[3];
        rand0_360(angles);
        //std::cout << " ------------------------------------------------ " << std::endl;
        //std::cout << "        Starting new particle orientation         " << std::endl;
        //std::cout << " ------------------------------------------------ " << std::endl;
        //std::cout << "Rotation angles: " << angles[0] << " " << angles[1] << " " << angles[2] << std::endl;
        //orientations_log << angles[0] << " " << angles[1] << " " << angles[2] << std::endl;
        crystal.rotateCrystal(angles[0], angles[1], angles[2], center);
        crystal.reset_structure();
        double random_magnetization_vector[3];
        marsaglia(random_magnetization_vector);
        for(int atom = 0; atom<(int)crystal.atoms.size(); atom++){
            if( crystal.atoms[atom].position == 0){
                crystal.atoms[atom].spinx = random_magnetization_vector[0];
                crystal.atoms[atom].spiny = random_magnetization_vector[1];
                crystal.atoms[atom].spinz = random_magnetization_vector[2];
            } else {
                crystal.atoms[atom].spinx = -random_magnetization_vector[0];
                crystal.atoms[atom].spiny = -random_magnetization_vector[1];
                crystal.atoms[atom].spinz = -random_magnetization_vector[2];
            }
        }
    //    Crystal crystal(structure_filename, dipoleInteractions, FeTT, FeOO, FeTO, FeOO_APB, anisotropyConstant, angles[0], angles[1], angles[2],macrocell_size, center);
        //crystal.update_temperature(TUpperLimit);
        int totalNumAtoms = (int)crystal.atoms.size();
  
        if(ZFC){
            
            // zero (small guide field) field cooling 
            //std::cout << "start cooling" << std::endl;
//            for(int j=0; j<temperature_arr_down.size();j++){
//                //crystal.update_temperature(temperature_arr_down[j]);
//                for(int k=0; k<(int)(steps*totalNumAtoms); k++){
//                    crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(cooling_field, temperature_arr_down[j]);
//                }
//              //  std::cout << temperature_arr_down[j] << std::endl;
//            }

//            for(int a=0; a<(int)(3000*totalNumAtoms); a++){
//                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(0.0, 0.01);
//            }
            //std::cout << "start ZFC measurement" << std::endl;
            for(int j=0; j<temperature_arr_up.size(); j++){
            
                // relaxation steps
                //crystal.update_temperature(temperature_arr_up[j]);
                for(int k = 0; k<(int)(steps*totalNumAtoms); k++){
                    crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field, temperature_arr_up[j]);
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
                        crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field, temperature_arr_up[j]);
                    
                        //omp_set_dynamic(0);
                        //# pragma omp parallel for reduction(+:mx_temp, my_temp, mz_temp) num_threads(3)
                        for(int a=0; a<(int)crystal.atoms.size(); a++){
                            mx_temp += MAGFE3*crystal.atoms[a].spinx;
                            my_temp += MAGFE3*crystal.atoms[a].spiny;
                            mz_temp += MAGFE3*crystal.atoms[a].spinz;
                        }
                    }
                    mxZFC[j] += (mx_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    myZFC[j] += (my_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    mzZFC[j] += (mz_temp/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                }
 
               // std::cout << temperature_arr_up[j] << std::endl;
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
//            for(int atom = 0; atom<(int)crystal.atoms.size(); atom++){
//                if( crystal.atoms[atom].position == 0){
//                    crystal.atoms[atom].spinx = 1.0;
//                    crystal.atoms[atom].spiny = 0.0;
//                    crystal.atoms[atom].spinz = 0.0;
//                } else {
//                    crystal.atoms[atom].spinx = -1.0;
//                    crystal.atoms[atom].spiny = 0.0;
//                    crystal.atoms[atom].spinz = 0.0;
//                }
//            }
            //std::cout << "start FC measurement" << std::endl;
            //crystal.reset_structure();
            for(int j=0; j<temperature_arr_down.size(); j++){
            //for(int j=0; j<temperature_arr_up.size(); j++){
                //crystal.update_temperature(temperature_arr_down[j]);
                
                // relaxation steps
                for(int k = 0; k< (int)(steps*totalNumAtoms); k++){
                    crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field, temperature_arr_down[j]);
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
                        crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(measurement_field, temperature_arr_down[j]);
                    
                        //omp_set_dynamic(0);
                        //# pragma omp parallel for reduction(+:mx_temp, my_temp, mz_temp) num_threads(3)
                        for(int a=0; a<(int)crystal.atoms.size(); a++){
                            mx_temp_fc += MAGFE3*crystal.atoms[a].spinx;
                            my_temp_fc += MAGFE3*crystal.atoms[a].spiny;
                            mz_temp_fc += MAGFE3*crystal.atoms[a].spinz;
                        }
                    }
                    mxFC[j] += (mx_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    myFC[j] += (my_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                    mzFC[j] += (mz_temp_fc/(double)(averaging_steps*totalNumAtoms))/(double)numOrientations;
                }

              //  std::cout << temperature_arr_down[j] << std::endl;
            } // end of measurement
           
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
            for(int currentStep = 0; currentStep < temperature_arr_down.size(); currentStep++){
                tSweepFc << temperature_arr_down[currentStep] << " " << mxFC[currentStep] << " " << myFC[currentStep] << " " << mzFC[currentStep] << "\n";
            }
        } // end of FC
        
    //std::cout << "orientation " << i+1 <<  " finished" << std::endl;
    }
}

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
                                      anisotropyConstant, ZFC, FC, center,lattice_a, lattice_b, lattice_c, sigma);
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
                                 int coolingSteps, double macrocell_size):dipoleInteractions(dipoleInteractions), steps(steps), numOrientations(numOrientations),   temperature(temperature),BUpperLimit(BUpperLimit), BLowerLimit(BLowerLimit), BstepSize(BstepSize), startTemp(startTemp), tempStep(tempStep), coolingField(coolingField), coolingSteps(coolingSteps), macrocell_size(macrocell_size) {
    
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

void MvsBMeasurement::fieldSweep(std::string output_dir, std::string structure_filename, double FeTT,double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant, double lattice_a,double lattice_b,double lattice_c,double center,double sigma){
    
    std::string output_filename = structure_filename.substr(structure_filename.find_last_of("/")+1);
    output_filename += "_Hysteresis_sim_"+std::to_string((int)(temperature));
    output_filename += "K_"+std::to_string((int)(steps))+"steps";
    output_filename += "_"+std::to_string((int)(numOrientations))+"or";
    output_filename += "_"+std::to_string((int)(coolingField))+"T";
    output_filename += "_"+std::to_string((int)(coolingSteps))+"cs";
    output_filename += "_sT"+std::to_string((int)(startTemp))+"K";
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
        
        bar.update(i);
        bar.set_status_text("orientation " + std::to_string(i+1));
        
        double angles[3];
        rand0_360(angles);
        //std::cout << "Rotation angles: " << angles[0] << " " << angles[1] << " " << angles[2] << std::endl;
        Crystal crystal(structure_filename, dipoleInteractions, FeTT, FeOO, FeTO, FeOO_APB, anisotropyConstant, angles[0], angles[1], angles[2],macrocell_size, center,lattice_a, lattice_b, lattice_c,sigma);
        
        int totalNumAtoms = int(crystal.atoms.size());
   
        //field cooling
        auto temperature_arr = arange<double>(startTemp,temperature,tempStep);
        
        for(int st=0; st<3000*totalNumAtoms; st++){
            crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(0.0, temperature_arr[0]);
        }
        //std::cout << "start cooling" << std::endl;
        for(int i=0; i<temperature_arr.size();i++){
            //crystal.update_temperature(temperature_arr[i]);
            //std::cout << "sigma: " << crystal.atoms[0].sigma << std::endl;
            for(int j=0; j<coolingSteps * totalNumAtoms; j++){
                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(coolingField, temperature_arr[i]);
            }
            //std::cout << temperature_arr[i] << std::endl;
        }
        //std::cout << "start field sweep" << std::endl;
        crystal.update_temperature(temperature);
        for(int j=0; j<=BnumberOfSteps; j++){
            
            for(int k = 0; k<steps*totalNumAtoms; k++){
                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(magneticField[j], temperature);
            }
            
            for(unsigned int l=0; l<=crystal.atoms.size(); l++){
                mX[j] += MAGFE3*crystal.atoms[l].spinx;
                mY[j] += MAGFE3*crystal.atoms[l].spiny;
                mZ[j] += MAGFE3*crystal.atoms[l].spinz;
            }
            //std::cout << magneticField[j] << std::endl;
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
        for(int currentStep = 0; currentStep <= BnumberOfSteps; currentStep++){
            bSweep << magneticField[currentStep] << " " << mX[currentStep] << " " << mY[currentStep] << " " << mZ[currentStep] << "\n";
        }
        
        //std::cout << "orientation " << i+1 <<  " finished" << std::endl;
    }
}

void run_MvsB(std::string output_dir, std::string structure_filename,std::string dipoleInteractions, int steps, int numOrientations, double temperature, double BUpperLimit, double BLowerLimit, double BstepSize, double coolingField, int coolingSteps, double startTemp, double tempStep, double FeTT, double FeOO, double FeTO, double FeOO_APB,double anisotropyConstant, double macrocell_size,double lattice_a, double lattice_b, double lattice_c, double center, double sigma){
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
                               anisotropyConstant, lattice_a,  lattice_b,  lattice_c,center, sigma);
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

void spinStructure::spinStructureMeasurement(double FeTT,double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant, double alpha, double beta, double gamma, double center, double lattice_a, double lattice_b, double lattice_c, double sigma){
    
    
    std::cout << "dev branch" << std::endl;
    Crystal crystal(structure_filename, dipoleInteractions,
                    FeTT, FeOO, FeTO, FeOO_APB, anisotropyConstant,
                    alpha, beta, gamma, macrocell_size, center,
                    lattice_a, lattice_b, lattice_c, sigma);
    crystal.update_temperature(temperature);

    int totalNumAtoms = int(crystal.atoms.size());

//    std::ofstream relaxation_oct;
//    relaxation_oct.open (filename+"relaxation_oct.dat", std::fstream::out);
//    std::ofstream relaxation_tet;
//    relaxation_tet.open (filename+"relaxation_tet.dat", std::fstream::out);

//    double magnetization_x_oct = 0.0;
//    double magnetization_y_oct = 0.0;
//    double magnetization_z_oct = 0.0;
//    double magnetization_x_tet = 0.0;
//    double magnetization_y_tet = 0.0;
//    double magnetization_z_tet = 0.0;
//    double magnetization_x = 0.0;
//    int numOct = 0;
//    int numTet = 0;
    
    ProgressBar bar;
    bar.set_bar_width(50);
    bar.fill_bar_progress_with("■");
    bar.fill_bar_remainder_with(" ");
    bar.update(0);

    //------Perform Monte Carlo Relaxation-----------------------------------------
    {
        for(int i =0; i < steps*totalNumAtoms; i++){
            
            
           // for(int k = 0; k < 1000000; k++){
                crystal.atoms[rand0_crystalAtoms(totalNumAtoms)].MonteCarloStep(magneticField, temperature);
                if(i % 100000 == 0){
                    bar.update((double)(i)/(double)(steps*totalNumAtoms));
                    bar.set_status_text("trial move " + std::to_string(i));
                }
        }
//            magnetization_x_oct = 0.0;
//            magnetization_y_oct = 0.0;
//            magnetization_z_oct = 0.0;
//            magnetization_x_tet = 0.0;
//            magnetization_y_tet = 0.0;
//            magnetization_z_tet = 0.0;
//            magnetization_x = 0.0;
//            numOct = 0;
//            numTet = 0;
//
//            for(int l=0; l<crystal.atoms.size(); l++){
//                magnetization_x += crystal.atoms[l].spinx;
//                if(crystal.atoms[l].position == 1){
//                    magnetization_x_oct += crystal.atoms[l].spinx;
//                    magnetization_y_oct += crystal.atoms[l].spiny;
//                    magnetization_z_oct += crystal.atoms[l].spinz;
//                    numOct++;
//                }
//                else{
//                    magnetization_x_tet += crystal.atoms[l].spinx;
//                    magnetization_y_tet += crystal.atoms[l].spiny;
//                    magnetization_z_tet += crystal.atoms[l].spinz;
//                    numTet++;
//                }
//            }
//
//            std::cout << magnetization_x << std::endl;
//
//            relaxation_oct << magnetization_x_oct/numOct << " " << magnetization_y_oct/numOct << " " << magnetization_z_oct/numOct << std::endl;
//            relaxation_tet << magnetization_x_tet/numTet << " " << magnetization_y_tet/numTet << " " << magnetization_z_tet/numTet << std::endl;
          
            //std::cout << i << std::endl;
            
        
        
        std::string output_filename = output_dir + structure_filename.substr(structure_filename.find_last_of("/")+1);
        output_filename += "_spin_structure_"+std::to_string((int)(temperature));
        output_filename += "K_"+std::to_string((int)(steps))+"MCS";
        output_filename += "_"+std::to_string((int)(magneticField))+"T";
        output_filename += "_dip"+dipoleInteractions;
        if(dipoleInteractions == "macrocell_method"){
//            std::string num_text = std::to_string(macrocell_size);
//            std::string rounded = num_text.substr(0, num_text.find(".")+3);
            
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
        structure << "# Orientation: " << alpha << ", " << beta << ", " << gamma << std::endl;
        structure << "# " << std::endl;
        structure << "# x      y        z      spinx    spiny    spinz  pos  APB" << std::endl;
        
        for(int i=0; i< crystal.atoms.size(); i++){
            structure << crystal.atoms[i].x << ", " << crystal.atoms[i].y << ", " << crystal.atoms[i].z <<", "<< crystal.atoms[i].spinx << ", " << crystal.atoms[i].spiny << ", " << crystal.atoms[i].spinz  << ", " << crystal.atoms[i].position << ", " << crystal.atoms[i].isApbAtom << "\n";
        }
    }
}

void run_spinstructure(std::string dipoleInteractions,int steps, double magneticField, double temperature, std::string output_dir, std::string structure_filename,double FeTT,double FeOO, double FeTO, double FeOO_APB,double anisotropyConstant, double alpha,double beta, double gamma, double macrocell_size, double center, double lattice_a, double lattice_b, double lattice_c, double sigma){
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
