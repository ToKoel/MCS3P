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
#include <sstream>
#include <functional>
#include <map>
#include <regex>
#include <algorithm>
#include <stdexcept>
#include "HelperStructs.hpp"

size_t split(const std::string, std::vector<std::string>, char);

class StructureFileParser{
public:
    static utility::StructureProperties parseStructureFile(std::string filename);
};

class CommandLineParser{
public:
    static utility::MeasurementSettings parseCommandline(char *argv[]);
};

#endif /* Parsers_hpp */
