//
//  Parsers.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.03.22.
//

#include "Parsers.hpp"

StructureProperties StructureFileParser::parseStructureFile(){
    std::ifstream inFile;
    inFile.open(filename);
    std::string line;
    int num = 0;
    while (std::getline(inFile, line)){
            ++num;
    }
    inFile.close();
    
    StructureProperties structureProperties;
    structureProperties.numberOfAtoms = num;
    structureProperties.positionVectors.resize(num);
    structureProperties.positionIDs.resize(num);
    structureProperties.isAPB.resize(num);

    inFile.open(filename);
    for(int i =0; i<num; i++){
        double xr,yr,zr;
        double posr,apbr,uc_xr,uc_yr,uc_zr;
        inFile >> xr >> yr >> zr >> posr >> apbr >> uc_xr >> uc_yr >> uc_zr ;
        structureProperties.positionVectors[i] = {xr, yr, zr};
        if(posr == 1){
            structureProperties.positionIDs[i] = StructuralPositions::kTetrahedral;
        } else {
            structureProperties.positionIDs[i] = StructuralPositions::kOctahedral;
        }
        if(apbr == 0){
            structureProperties.isAPB[i] = false;
        } else {
            structureProperties.isAPB[i] = true;
        }
    }
    inFile.close();
    
    return structureProperties;
}
