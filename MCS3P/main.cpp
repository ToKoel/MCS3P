#include "AtomTests.hpp"
#include "EndToEndTests.hpp"
#include "Measurement.hpp"
#include "ParserTests.hpp"
#include "randNumGenerator.hpp"
#include "testing.hpp"
#include "gtest/gtest.h"
#include <iostream>

int main(int argc, char *argv[]) {
  seed250(SEED);

  if (argc > 1) {
      utility::MeasurementSettings measurementSettings =
        CommandLineParser::parseCommandline(argv);
    Measurement measurement(measurementSettings);

    switch (measurementSettings.measurementType) {
        case utility::MeasurementType::kNone:
      std::cout << "no measurement selected\n";
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
