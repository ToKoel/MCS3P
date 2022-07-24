#ifndef Measurement_hpp
#define Measurement_hpp

#include <vector>
#include <string>
#include "ProgressBar.hpp"
#include "Crystal.hpp"
#include "HelperStructs.hpp"
#include "FileIO.hpp"


class Measurement{
public:
    utility::MeasurementSettings measurementSettings;
    fileIO::FileHandler fileHandler{measurementSettings};
    Measurement(utility::MeasurementSettings measurementSettings): measurementSettings(measurementSettings) {};
    void run_MvsT();
    void run_MvsB();
    void run_spinstructure();
    
    void temperatureSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<utility::LinalgVector>& magnetizationArray);
    void fieldSweep(Crystal& crystal, std::vector<double>& temperatureArray, std::vector<double>& fieldArray, std::vector<utility::LinalgVector>& magnetizationArray);
};

#endif /* Measurement_hpp */
