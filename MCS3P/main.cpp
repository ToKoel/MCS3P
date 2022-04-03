#include "randNumGenerator.hpp"
#include "Measurement.hpp"
#include "testing.hpp"

int main(int argc, char* argv[]) {
    seed250(SEEED);
    
    MeasurementSettings measurementSettings = CommandLineParser::parseCommandline(argv[1]);
    Measurement measurement(measurementSettings);
 
    switch(measurementSettings.measurementType){
        case MeasurementType::kNone:
            break;
        case MeasurementType::kSpinStructure:
            measurement.run_spinstructure();
            break;
        case MeasurementType::kMvsH:
            measurement.run_MvsB();
            break;
        case MeasurementType::kMvsT:
            measurement.run_MvsT();
            break;
        case MeasurementType::kTest:
            break;
    }
    
    return 0;
}
