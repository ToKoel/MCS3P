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
template<typename T>
std::string to_string(const T& obj, int prec)
{
    std::ostringstream stream;
    stream.precision(prec);
    stream<<std::fixed;
    stream << obj;
    std::string stringVar = stream.str();
    std::replace(stringVar.begin(), stringVar.end(), '.', 'd');
    return stringVar;
}

std::string Measurement::generateOutputHeader(){
    std::stringstream stream;
    switch(measurementSettings.measurementType){
        case MeasurementType::kMvsT:
            stream << "# structure file: " << measurementSettings.structurePath << "\n";
            stream << "# dipole interactions: " << DipoleInteractionTypes[static_cast<int>(measurementSettings.dipoleInteractionHandling)] << "\n";
            stream << "# steps: " << measurementSettings.monteCarloSteps << "\n";
            stream << "# number of orientations: " << measurementSettings.numOrientations  << "\n";
            stream << "# measurement field: " << measurementSettings.measurementField << "\n";
            stream << "# cooling field: " << measurementSettings.coolingField << "\n";
            stream << "# upper temperature limit: " << measurementSettings.upperTemperatureLimit << "\n";
            stream << "# lower temperature limit: " << measurementSettings.lowerTemperatureLimit << "\n";
            stream << "# temperature step size: " << measurementSettings.temperatureStepSize << "\n";
            stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
            stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
            stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
            stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB << "\n";
            stream << "# lattice parameters: " << measurementSettings.latticeParameters.a;
            stream << ", " << measurementSettings.latticeParameters.b;
            stream << ", " << measurementSettings.latticeParameters.c << "\n";
            stream << "# sigma: " << measurementSettings.sigma << "\n";
            stream << "# " << "\n";
            stream << "# field (T)      M_x(B)      M_y(B)      M_z(B)" << "\n";
            break;
            
        case MeasurementType::kMvsH:
            stream << "# structure file: " << measurementSettings.structurePath << "\n";
            stream << "# dipole interactions: " << DipoleInteractionTypes[static_cast<int>(measurementSettings.dipoleInteractionHandling)] << "\n";
            stream << "# steps: " << measurementSettings.monteCarloSteps << "\n";
            stream << "# number of orientations: " << measurementSettings.numOrientations  << "\n";
            stream << "# temperature: " << measurementSettings.temperature << "\n";
            stream << "# B_upper: " << measurementSettings.upperFieldLimit << "\n";
            stream << "# B_lower: " << measurementSettings.lowerFieldLimit << "\n";
            stream << "# B_step: " << measurementSettings.fieldStepSize << "\n";
            stream << "# cooling_field: " << measurementSettings.coolingField << "\n";
            stream << "# start_temperature: " << measurementSettings.startingTemperature << "\n";
            stream << "# temperature_step: " << measurementSettings.temperatureStepSize << "\n";
            stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
            stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
            stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
            stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB << "\n";
            stream << "# lattice parameters: " << measurementSettings.latticeParameters.a;
            stream << ", " << measurementSettings.latticeParameters.b;
            stream << ", " << measurementSettings.latticeParameters.c << "\n";
            //stream << "# total number of iron atoms: " << totalNumAtoms << std::endl;
            stream << "# " << std::endl;
            stream << "# field (T)      M_x(B)      M_y(B)      M_z(B)" << "\n";
            break;
            
        case MeasurementType::kSpinStructure:
            stream << "# structure_file: " << measurementSettings.structurePath << "\n";
            stream << "# dipole_interactions: " << DipoleInteractionTypes[static_cast<int>(measurementSettings.dipoleInteractionHandling)] << "\n";
            if(measurementSettings.dipoleInteractionHandling == DipoleInteractions::kMacrocellMethod){
                stream << "# macrocell_size: " << measurementSettings.macrocellSize << "\n";
            }
            stream << "# steps: " << measurementSettings.monteCarloSteps << std::endl;
            stream << "# Number of atoms: " << measurementSettings.totalNumAtoms << "\n";
            stream << "# temperature: " << measurementSettings.temperature << "\n";
            stream << "# field: " << measurementSettings.measurementField << "\n";
            stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
            stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
            stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
            stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB << "\n";
            stream << "# anisotropy constant: " << measurementSettings.anisotropyConstant << "\n";
            stream << "# lattice parameters: " << measurementSettings.latticeParameters.a;
            stream << ", " << measurementSettings.latticeParameters.b;
            stream << ", " << measurementSettings.latticeParameters.c << "\n";
            stream << "# Orientation: " << measurementSettings.angles.alpha << ", ";
            stream << measurementSettings.angles.beta << ", " << measurementSettings.angles.gamma << "\n";
            stream << "# " << "\n";
            stream << "# x      y        z      spinx    spiny    spinz  pos  APB" << "\n";
            break;
            
        case MeasurementType::kNone:
            break;
            
        case MeasurementType::kTest:
            break;
    }
    
    return stream.str();
}

std::string Measurement::generateOutputFilename(){
    std::string outputFilename = measurementSettings.outputPath;
    outputFilename += measurementSettings.structurePath.substr(measurementSettings.structurePath.find_last_of("/")+1);
    switch(measurementSettings.measurementType){
        case MeasurementType::kMvsT:
            outputFilename += "_MvsT_sim_cF" + to_string(measurementSettings.coolingField, 2);
            outputFilename += "T_mF" + to_string(measurementSettings.measurementField, 3);
            outputFilename += "T_" + to_string(measurementSettings.monteCarloSteps, 1) + "steps";
            outputFilename += "_" + to_string(measurementSettings.averagingSteps, 1) + "avsteps";
            outputFilename += "_" + to_string(measurementSettings.numOrientations, 1) + "or";
            outputFilename += "_" + to_string(measurementSettings.sigma, 2) + "sig";
            break;
            
        case MeasurementType::kMvsH:
            outputFilename += "_Hysteresis_sim_" + to_string(measurementSettings.temperature, 1);
            outputFilename += "K_" + to_string(measurementSettings.monteCarloSteps, 1) + "steps";
            outputFilename += "_" + to_string(measurementSettings.numOrientations, 1) + "or";
            outputFilename += "_" + to_string(measurementSettings.coolingField, 2) + "T";
            outputFilename += "_" + to_string(measurementSettings.coolingSteps, 2) + "cs";
            outputFilename += "_sT"+ to_string(measurementSettings.startingTemperature, 1) + "K";
            outputFilename += "_dip"+ DipoleInteractionTypes[static_cast<int>(measurementSettings.dipoleInteractionHandling)];
            break;
            
        case MeasurementType::kSpinStructure:
            outputFilename += "_spin_structure_"+ to_string(measurementSettings.temperature, 1);
            outputFilename += "K_" + to_string(measurementSettings.monteCarloSteps, 1) +"MCS";
            outputFilename += "_"+ to_string(measurementSettings.measurementField, 1) +"T";
            outputFilename += "_dip"+ DipoleInteractionTypes[static_cast<int>(measurementSettings.dipoleInteractionHandling)];
            break;
            
        case MeasurementType::kNone:
            break;
            
        case MeasurementType::kTest:
            break;
    }
    
    if(measurementSettings.dipoleInteractionHandling == DipoleInteractions::kMacrocellMethod){
        outputFilename += "_" + to_string(measurementSettings.macrocellSize, 4) + "mcsize";
    }
    
    return outputFilename;
}


void Measurement::temperatureSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<LinalgVector>& magnetizationArray){
    auto totalNumAtoms = crystal.atoms.size();
    for(auto i=0; i<temperatureArray.size(); i++){
        crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps,
                                       measurementSettings.measurementField,
                                       temperatureArray[i]);
    
        crystal.performMonteCarloSteps(measurementSettings.averagingSteps,
                                       measurementSettings.measurementField,
                                       temperatureArray[i]);

        magnetizationArray[i] += crystal.getNetSpinVector()/(static_cast<double>(totalNumAtoms*measurementSettings.numOrientations));
    }
}

void Measurement::fieldSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<double>& fieldArray, std::vector<LinalgVector>& magnetizationArray){
    auto totalNumAtoms = crystal.atoms.size();
    
    for(int i=0; i<temperatureArray.size();i++){
        crystal.performMonteCarloSteps(measurementSettings.coolingSteps * totalNumAtoms,
                                       measurementSettings.coolingField,
                                       temperatureArray[i]);
    }
    
    for(int i=0; i<fieldArray.size();i++){
        crystal.performMonteCarloSteps(measurementSettings.coolingSteps * totalNumAtoms,
                                       fieldArray[i],
                                       measurementSettings.temperature);
        
        magnetizationArray[i] += crystal.getNetSpinVector()/(static_cast<double>(totalNumAtoms*measurementSettings.numOrientations));
    }
}

void Measurement::writeSweepResultsToFile(MeasurementMode mode, std::vector<double>& variableArray, std::vector<LinalgVector>& magnetizationArray){
    std::ofstream sweep;
    switch(mode){
        case MeasurementMode::ZFC:
            sweep.open(generateOutputFilename() + "_ZFC.txt", std::fstream::out);
            break;
        case MeasurementMode::FC:
            sweep.open(generateOutputFilename() + "_FC.txt", std::fstream::out);
            break;
        case MeasurementMode::Field: 
            sweep.open(generateOutputFilename() + ".txt", std::fstream::out);
            break;
    }
    sweep << generateOutputHeader();
    for(int currentStep = 0; currentStep < variableArray.size(); currentStep++){
        sweep << variableArray[currentStep] << " " << magnetizationArray[currentStep].x << " " << magnetizationArray[currentStep].y << " " << magnetizationArray[currentStep].z << "\n";
    }
    sweep.close();
}

void Measurement::writeSpinStructureResultsToFile(Crystal& crystal){
    std::ofstream structure;
    structure.open(generateOutputFilename() + ".txt", std::fstream::out);
    structure << generateOutputHeader();
    for(auto atom: crystal.atoms){
        structure << atom.positionVector.x << ", "
                  << atom.positionVector.y << ", "
                  << atom.positionVector.z << ", "
                  << atom.spinVector.x << ", "
                  << atom.spinVector.y << ", "
                  << atom.spinVector.z << ", "
                  << StructurePositionTypes[static_cast<int>(atom.structuralPositionID)] << ", "
                  << atom.isApbAtom << "\n";
    }
    structure.close();
}

void Measurement::run_MvsT(){
    std::string outputFileName = generateOutputFilename();
    std::cout << "\nOutput filename: " << outputFileName << "\n";
    std::cout << "\n ------- starting M vs T measurement -------\n" << "\n";
    
    auto temperatureStepsCooling = arange<double>(measurementSettings.upperTemperatureLimit,
                                               measurementSettings.lowerTemperatureLimit,
                                               measurementSettings.temperatureStepSize);
    auto temperatureStepsHeating = arange<double>(measurementSettings.lowerTemperatureLimit,
                                             measurementSettings.upperTemperatureLimit,
                                             measurementSettings.temperatureStepSize);
    
    std::vector<LinalgVector> magnetizationZFC;
    std::vector<LinalgVector> magnetizationFC;
    
    magnetizationZFC.resize(temperatureStepsHeating.size());
    magnetizationFC.resize(temperatureStepsCooling.size());
    
    ProgressBar bar;
    
    Crystal crystal(measurementSettings);

    for(int i = 0; i < measurementSettings.numOrientations; i++){
        bar.step(static_cast<double>(i)/measurementSettings.numOrientations, "orientation " + std::to_string(i+1));
        
        crystal.randomOrientation();
        crystal.alignAlongRandomVector();

        if(measurementSettings.ZFC){
            temperatureSweep(crystal, temperatureStepsHeating, magnetizationZFC);
            writeSweepResultsToFile(MeasurementMode::ZFC, temperatureStepsHeating, magnetizationZFC);
        }
        
        if(measurementSettings.FC){
            temperatureSweep(crystal, temperatureStepsCooling, magnetizationFC);
            writeSweepResultsToFile(MeasurementMode::FC, temperatureStepsCooling, magnetizationFC);
        }
    }
}


void Measurement::run_MvsB(){
    std::string outputFileName = generateOutputFilename();
    std::cout << "\nOutput filename: " << outputFileName << "\n";
    std::cout << "\n ------- starting M vs B measurement -------\n" << "\n";
    
    std::vector<double> fieldArray       = arange<double>(0,
                                                          measurementSettings.upperFieldLimit,
                                                          measurementSettings.fieldStepSize);
    std::vector<double> fieldArrayUp     = arange<double>(measurementSettings.lowerFieldLimit,
                                                          measurementSettings.upperFieldLimit,
                                                          measurementSettings.fieldStepSize);
    std::vector<double> fieldArrayDown   = arange<double>(measurementSettings.upperFieldLimit,
                                                          measurementSettings.lowerFieldLimit,
                                                          measurementSettings.fieldStepSize);
    std::vector<double> temperatureArray = arange<double>(measurementSettings.startingTemperature,
                                                          measurementSettings.temperature,
                                                          measurementSettings.temperatureStepSize);
    fieldArray.insert( fieldArray.end(), fieldArrayDown.begin(), fieldArrayDown.end() );
    fieldArray.insert( fieldArray.end(), fieldArrayUp.begin(), fieldArrayUp.end() );
    fieldArray.push_back(measurementSettings.upperFieldLimit);
    std::vector<LinalgVector> magnetization;
    magnetization.resize(fieldArray.size());
    
    ProgressBar bar;
    for(int i = 0; i < measurementSettings.numOrientations; i++){
        bar.step(static_cast<double>(i)/measurementSettings.numOrientations, "orientation " + std::to_string(i+1));
        Crystal crystal(measurementSettings);
        crystal.randomOrientation();
        fieldSweep(crystal, temperatureArray, fieldArray, magnetization);
    }
    writeSweepResultsToFile(MeasurementMode::Field, fieldArray, magnetization);
}



void Measurement::run_spinstructure(){
    std::string outputFileName = generateOutputFilename();
    std::cout << "\nOutput filename: " << outputFileName << "\n";
    std::cout << "\n ------- starting spin structure measurement -------\n" << "\n";
    
    Crystal crystal(measurementSettings);

    //simulated annealing
    std::vector<double> coolingArray = arange(measurementSettings.startingTemperature,
                                              measurementSettings.temperature,
                                              measurementSettings.coolingSteps);
    coolingArray.push_back(measurementSettings.temperature);
    for(double temperatureStep: coolingArray){
        crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps/100,
                                       measurementSettings.measurementField,
                                       temperatureStep);
    }
    
    // relaxation
    crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps,
                                   measurementSettings.measurementField,
                                   measurementSettings.temperature);
    
    writeSpinStructureResultsToFile(crystal);
}
