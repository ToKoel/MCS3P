//
//  EndToEndTests.h
//  MCS3P
//
//  Created by Tobias Köhler on 17.04.22.
//

#ifndef EndToEndTests_h
#define EndToEndTests_h

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <filesystem>
#include <fstream>
#include <string>
#include "Parsers.hpp"
#include "Crystal.hpp"




//TEST(EndToEndTest, atomicSpinsCanBeAligned){
//    std::filesystem::path dir (std::filesystem::temp_directory_path());
//    std::filesystem::path structureFile ("D8_structure_noAPB");
//    std::filesystem::path structurePath = dir / structureFile;
//    std::filesystem::path origin = std::filesystem::current_path();
//    std::filesystem::path originStructureFile ("TestFiles/D8_structure_noAPB");
//    std::filesystem::path originStructurePath = origin / originStructureFile;
//    
//    bool worked = copy_file( originStructurePath, structurePath);
//    std::cout << worked << std::endl;
//    
//    utility::MeasurementSettings measurementSettings;
//    measurementSettings.structurePath = structurePath.string();
//    Crystal crystal(measurementSettings);
//    double center = 4.5;
//    
//    std::cout << crystal.atoms[0].positionVector.y << std::endl;
//    
////    for(auto& atom : crystal.atoms){
////        if (atom.positionVector.y > center){
////            if (atom.structuralPositionID == StructuralPositions::kOctahedral){
////                atom.spinVector = {-1.0, 0.0, 0.0};
////            } else{
////                atom.spinVector = {1.0, 0.0, 0.0};
////            }
////        } else if (atom.positionVector.y <= center){
////            if (atom.structuralPositionID == StructuralPositions::kOctahedral){
////                atom.spinVector = {1.0, 0.0, 0.0};
////            } else{
////                atom.spinVector = {-1.0, 0.0, 0.0};
////            }
////        }
////    }
////
////    ASSERT_DOUBLE_EQ(crystal.atoms[1].spinVector.x, -1.0);
//}


TEST(EndToEndTest, parserGeneratesMeasurementSettings){
    std::string file = "../TestFiles/testSettings.json";
    char* testFile[] = {file.data()};
    utility::MeasurementSettings settings = CommandLineParser::parseCommandline(testFile);
    
    ASSERT_EQ(settings.monteCarloSteps, 5);
    ASSERT_EQ(settings.particleSize, 9);
    ASSERT_DOUBLE_EQ(settings.temperature, 0.004);
}

#endif /* EndToEndTests_h */
