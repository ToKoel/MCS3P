//
//  Parsers.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.03.22.
//

#include "Parsers.hpp"

utility::StructureProperties
StructureFileParser::parseStructureFile(std::string filename) {
  std::ifstream inFile;
  inFile.open(filename);
  std::string line;
  size_t num = 0;
  while (std::getline(inFile, line)) {
    ++num;
  }
  inFile.close();

    utility::StructureProperties structureProperties;
  structureProperties.numberOfAtoms = num;
  structureProperties.positionVectors.resize(num);
  structureProperties.positionIDs.resize(num);
  structureProperties.isAPB.resize(num);

  inFile.open(filename);
  for (size_t i = 0; i < static_cast<size_t>(num); i++) {
    double xr, yr, zr;
    double posr, apbr, uc_xr, uc_yr, uc_zr;
    inFile >> xr >> yr >> zr >> posr >> apbr >> uc_xr >> uc_yr >> uc_zr;
    structureProperties.positionVectors[i] = {xr, yr, zr};
    if (posr == 1) {
        structureProperties.positionIDs[i] = utility::StructuralPositions::kTetrahedral;
    } else {
        structureProperties.positionIDs[i] = utility::StructuralPositions::kOctahedral;
    }
    if (apbr == 0) {
      structureProperties.isAPB[i] = false;
    } else {
      structureProperties.isAPB[i] = true;
    }
  }
  inFile.close();

  return structureProperties;
}

utility::MeasurementSettings CommandLineParser::parseCommandline(char *argv[]) {
  std::string measurementSettingsFile = argv[0];
  std::ifstream file(measurementSettingsFile);
  std::string str;
  std::string delimiter = ":";
  std::regex regex(".*:.*");

  if (!file.good()) {
    throw std::invalid_argument("file doesn't exist");
  }

    utility::MeasurementSettings measurementSettings = {};

  while (std::getline(file, str)) {
    if (std::regex_match(str, regex)) {
      str.erase(std::remove(str.begin(), str.end(), ','), str.end());
      std::string token = str.substr(5, str.find(delimiter) - 6);
      std::string value = str.substr(str.find(delimiter) + 2);
      if (token == "measurement") {
        if (value == "M vs B") {
            measurementSettings.measurementType = utility::MeasurementType::kMvsH;
        } else if (value == "M vs T") {
            measurementSettings.measurementType = utility::MeasurementType::kMvsT;
        } else if (value == "spin structure") {
            measurementSettings.measurementType = utility::MeasurementType::kSpinStructure;
        } else if (value == "testing") {
            measurementSettings.measurementType = utility::MeasurementType::kTest;
        } else {
            measurementSettings.measurementType = utility::MeasurementType::kNone;
        }
      } else if (token == "dipole") {
        if (value == "brute_force") {
          measurementSettings.dipoleInteractionHandling =
            utility::DipoleInteractions::kBruteForce;
        } else if (value == "macrocell") {
          measurementSettings.dipoleInteractionHandling =
            utility::DipoleInteractions::kMacrocellMethod;
        } else {
          measurementSettings.dipoleInteractionHandling =
            utility::DipoleInteractions::kNoInteractions;
        }
      } else if (token == "output") {
        measurementSettings.outputPath = value;
      } else if (token == "structure_path") {
        measurementSettings.structurePath = value;
      } else if (token == "num_particles") {
        measurementSettings.numParticles = std::stoi(value);
      } else if (token == "particle_size") {
        measurementSettings.particleSize = std::stoi(value);
      } else if (token == "meas_field") {
        measurementSettings.measurementField = std::stod(value);
      } else if (token == "temperature") {
        measurementSettings.temperature = std::stod(value);
      } else if (token == "steps") {
        measurementSettings.monteCarloSteps = std::stoi(value);
      } else if (token == "av_steps") {
        measurementSettings.averagingSteps = std::stoi(value);
      } else if (token == "num_or") {
        measurementSettings.numOrientations = std::stoi(value);
      } else if (token == "cool_field") {
        measurementSettings.coolingField = std::stod(value);
      } else if (token == "Tupper") {
        measurementSettings.upperTemperatureLimit = std::stod(value);
      } else if (token == "Tlower") {
        measurementSettings.lowerTemperatureLimit = std::stod(value);
      } else if (token == "T_step") {
        measurementSettings.temperatureStepSize = std::stod(value);
      } else if (token == "alpha") {
        measurementSettings.angles.alpha = std::stod(value);
      } else if (token == "beta") {
        measurementSettings.angles.beta = std::stod(value);
      } else if (token == "gamma") {
        measurementSettings.angles.gamma = std::stod(value);
      } else if (token == "macrocell_size") {
        measurementSettings.macrocellSize = std::stod(value);
      } else if (token == "anisotropy_constant") {
        measurementSettings.anisotropyConstant = std::stod(value);
      } else if (token == "ZFC") {
        if (value == "false") {
          measurementSettings.ZFC = false;
        } else {
          measurementSettings.ZFC = true;
        }
      } else if (token == "FC") {
        if (value == "false") {
          measurementSettings.FC = false;
        } else {
          measurementSettings.FC = true;
        }
      } else if (token == "APB") {
        if (value == "false") {
          measurementSettings.APB = false;
        } else {
          measurementSettings.APB = true;
        }
      } else if (token == "Bupper") {
        measurementSettings.upperFieldLimit = std::stod(value);
      } else if (token == "Blower") {
        measurementSettings.lowerFieldLimit = std::stod(value);
      } else if (token == "Bstep") {
        measurementSettings.fieldStepSize = std::stod(value);
      } else if (token == "start_temperature") {
        measurementSettings.startingTemperature = std::stod(value);
      } else if (token == "cooling_steps") {
        measurementSettings.coolingSteps = std::stod(value);
      } else if (token == "sigma") {
        measurementSettings.sigma = std::stod(value);
      } else if (token == "APB_constant") {
        measurementSettings.exchangeConstants.FeOO_APB = std::stod(value);
      } else if (token == "Fe_TT") {
        measurementSettings.exchangeConstants.FeTT = std::stod(value);
      } else if (token == "Fe_OO") {
        measurementSettings.exchangeConstants.FeOO = std::stod(value);
      } else if (token == "Fe_TO") {
        measurementSettings.exchangeConstants.FeTO = std::stod(value);
      } else if (token == "lattice_a") {
        measurementSettings.latticeParameters.a = std::stod(value);
      } else if (token == "lattice_b") {
        measurementSettings.latticeParameters.b = std::stod(value);
      } else if (token == "lattice_c") {
        measurementSettings.latticeParameters.c = std::stod(value);
      } else if (token == "nearest_neighbour_distance") {
        measurementSettings.nearestNeighbourDistance = std::stod(value);
      }
    }
  }
  return measurementSettings;
}
