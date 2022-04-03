//
//  Parsers.hpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.03.22.
//

#ifndef Parsers_hpp
#define Parsers_hpp

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <regex>
#include <algorithm>
#include "HelperStructs.hpp"

class StructureFileParser{
public:
    static StructureProperties parseStructureFile(std::string filename);
};

class CommandLineParser{
public:
    static MeasurementSettings parseCommandline(std::string argv);
};

#endif /* Parsers_hpp */
