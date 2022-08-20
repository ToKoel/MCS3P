#include "AtomTests.hpp"
#include "EndToEndTests.hpp"
#include "Measurement.hpp"
#include "ParserTests.hpp"
#include "randNumGenerator.hpp"
#include "testing.hpp"
#include "gtest/gtest.h"
#include <iostream>

int main(int argc, char* argv[]) {
  if (argc > 1) {
      utility::MeasurementSettings measurementSettings =
        CommandLineParser::parseCommandline(argv);
    seed250(measurementSettings.seed);
    Measurement measurement(measurementSettings);

    switch (measurementSettings.measurementType) {
        case utility::MeasurementType::kNone:
      break;
        case utility::MeasurementType::kSpinStructure:
      measurement.run_spinstructure();
      break;
        case utility::MeasurementType::kMvsH:
      measurement.run_MvsB();
      break;
        case utility::MeasurementType::kMvsT:
      measurement.run_MvsT();
      break;
        case utility::MeasurementType::kTest:
      break;
    }
  } else {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

  return 0;
}
