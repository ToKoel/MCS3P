//
//  FileIO.h
//  MCS3P
//
//  Created by Tobias Köhler on 23.07.22.
//

#ifndef FileIO_hpp
#define FileIO_hpp

#include <string>
#include <sstream>
#include <fstream>
#include "HelperStructs.hpp"
#include "Crystal.hpp"

namespace fileIO {

class FileHandler {
private:
    utility::MeasurementSettings measurementSettings;
public:
    FileHandler(utility::MeasurementSettings measurementSettings);
    std::string generateOutputHeader();
    std::string generateOutputFilename();
    void writeSweepResultsToFile(utility::MeasurementMode mode, std::vector<double> &variableArray,
                            std::vector<utility::LinalgVector> &magnetizationArray);
    void writeSpinStructureResultsToFile(Crystal &crystal);
};



}

#endif /* FileIO_hpp */
