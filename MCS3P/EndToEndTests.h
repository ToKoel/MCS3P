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
#include "Parsers.hpp"

TEST(EndToEndTest, parserGeneratesMeasurementSettings){
    std::string file = "/Users/tobiaskohler/Desktop/MCS3P/MCS3P/TestFiles/testSettings.json";
    char* testFile[] = {file.data()};
    MeasurementSettings settings = CommandLineParser::parseCommandline(testFile);
    
    ASSERT_EQ(settings.monteCarloSteps, 5);
    ASSERT_EQ(settings.particleSize, 9);
    ASSERT_DOUBLE_EQ(settings.temperature, 0.004);
}

#endif /* EndToEndTests_h */
