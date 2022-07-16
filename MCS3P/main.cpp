#include <iostream>
#include "randNumGenerator.hpp"
#include "Measurement.hpp"
#include "testing.hpp"
#include "gtest/gtest.h"
#include "ParserTests.hpp"
#include "AtomTests.hpp"


int main(int argc, char* argv[]) {
    seed250(SEED);
    
    if(argc > 1){
        MeasurementSettings measurementSettings = CommandLineParser::parseCommandline(argv);
        Measurement measurement(measurementSettings);
     
        switch(measurementSettings.measurementType){
            case MeasurementType::kNone:
                std::cout << "no measurement selected" << std::endl;
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
    }
    else{
        testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
    
    return 0;
}
