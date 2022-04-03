#ifndef Measurement_hpp
#define Measurement_hpp

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include "ProgressBar.hpp"
#include "Crystal.hpp"


enum class MeasurementMode{
    ZFC,
    FC,
    Field
};

class Measurement{
public:
    MeasurementSettings measurementSettings;
    Measurement(MeasurementSettings measurementSettings): measurementSettings(measurementSettings) {};
    void run_MvsT();
    void run_MvsB();
    void run_spinstructure();
    
    std::string generateOutputHeader();
    std::string generateOutputFilename();
    void temperatureSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<LinalgVector>& magnetizationArray);
    void fieldSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<double>& fieldArray, std::vector<LinalgVector>& magnetizationArray);
    void writeSweepResultsToFile(MeasurementMode mode, std::vector<double>& temperatureArray, std::vector<LinalgVector>& magnetizationArray);
    void writeSpinStructureResultsToFile(Crystal& crystal);
};

#endif /* Measurement_hpp */
