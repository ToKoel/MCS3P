////
////  ParserTest.hpp
////  MCS3P
////
////  Created by Tobias Köhler on 12.04.22.
////
//
//#ifndef ParserTest_hpp
//#define ParserTest_hpp
//
//#include <gtest/gtest.h>
//#include <gmock/gmock.h>
//#include "Parsers.hpp"
//
//TEST(ParserTest, parserGeneratesMeasurementSettings){
//    std::string file = "/Users/tobiaskohler/Desktop/MCS3P/MCS3P/TestFiles/testSettings.json";
//    char* testFile[] = {file.data()};
//    MeasurementSettings settings = CommandLineParser::parseCommandline(testFile);
//    
//    ASSERT_EQ(settings.monteCarloSteps, 5);
//    ASSERT_EQ(settings.particleSize, 9);
//    ASSERT_DOUBLE_EQ(settings.temperature, 0.004);
//}
//
//TEST(ParserTest, parserGeneratesStructure){
//    std::string file = "/Users/tobiaskohler/Desktop/MCS3P/MCS3P/TestFiles/testStructure.txt";
//    StructureProperties structure = StructureFileParser::parseStructureFile(file);
//    
//    ASSERT_EQ(structure.numberOfAtoms, 6);
//    ASSERT_FALSE(structure.isAPB[4]);
//    ASSERT_TRUE(structure.isAPB[5]);
//    ASSERT_EQ(structure.positionVectors[0].x, 0.75000);
//    ASSERT_EQ(structure.positionVectors[0].y, 3.00000);
//    ASSERT_EQ(structure.positionVectors[0].z, 3.12500);
//    ASSERT_THAT(structure.positionIDs[0], StructuralPositions::kTetrahedral);
//}
//
//#endif /* ParserTest_hpp */
