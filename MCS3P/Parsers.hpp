//
//  Parsers.hpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.03.22.
//

#ifndef Parsers_hpp
#define Parsers_hpp

#include <vector>
#include "HelperStructs.hpp"
#include <string>
#include <fstream>


struct StructureProperties{
    int numberOfAtoms;
    std::vector<LinalgVector> positionVectors;
    std::vector<StructuralPositions> positionIDs;
    std::vector<bool> isAPB;
};

class StructureFileParser{
    std::string filename;
public:
    StructureProperties parseStructureFile();
    
    StructureFileParser(std::string filename): filename(filename){};
};


#endif /* Parsers_hpp */
