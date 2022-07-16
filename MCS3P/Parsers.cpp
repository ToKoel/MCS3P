//
//  Parsers.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.03.22.
//

#include "Parsers.hpp"

StructureProperties
StructureFileParser::parseStructureFile(std::string filename) {
  std::ifstream inFile;
  inFile.open(filename);
  std::string line;
  int num = 0;
  while (std::getline(inFile, line)) {
    ++num;
  }
  inFile.close();

  StructureProperties structureProperties;
  structureProperties.numberOfAtoms = num;
  structureProperties.positionVectors.resize(num);
  structureProperties.positionIDs.resize(num);
  structureProperties.isAPB.resize(num);

  inFile.open(filename);
  for (int i = 0; i < num; i++) {
    double xr, yr, zr;
    double posr, apbr, uc_xr, uc_yr, uc_zr;
    inFile >> xr >> yr >> zr >> posr >> apbr >> uc_xr >> uc_yr >> uc_zr;
    structureProperties.positionVectors[i] = {xr, yr, zr};
    if (posr == 1) {
      structureProperties.positionIDs[i] = StructuralPositions::kTetrahedral;
    } else {
      structureProperties.positionIDs[i] = StructuralPositions::kOctahedral;
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

size_t split(const std::string &txt, std::vector<std::string> &strs, char ch) {
  size_t pos = txt.find(ch);
  size_t initialPos = 0;
  strs.clear();

  // Decompose statement
  while (pos != std::string::npos) {
    strs.push_back(txt.substr(initialPos, pos - initialPos));
    initialPos = pos + 1;

    pos = txt.find(ch, initialPos);
  }

  // Add the last one
  strs.push_back(
      txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

  return strs.size();
}

MeasurementSettings CommandLineParser::parseCommandline(char *argv[]) {
  std::string measurementSettingsFile = argv[0];
  std::ifstream file(measurementSettingsFile);
  std::string str;
  std::string delimiter = ":";
  std::regex regex(".*:.*");

  if (!file.good()) {
    throw std::invalid_argument("file doesn't exist");
  }

  MeasurementSettings measurementSettings = {};

  while (std::getline(file, str)) {
    if (std::regex_match(str, regex)) {
      str.erase(std::remove(str.begin(), str.end(), ','), str.end());
      std::string token = str.substr(5, str.find(delimiter) - 6);
      std::string value = str.substr(str.find(delimiter) + 2);
      if (token == "measurement") {
        if (value == "M vs B") {
          measurementSettings.measurementType = MeasurementType::kMvsH;
        } else if (value == "M vs T") {
          measurementSettings.measurementType = MeasurementType::kMvsT;
        } else if (value == "spin structure") {
          measurementSettings.measurementType = MeasurementType::kSpinStructure;
        } else if (value == "testing") {
          measurementSettings.measurementType = MeasurementType::kTest;
        } else {
          measurementSettings.measurementType = MeasurementType::kNone;
        }
      } else if (token == "dipole") {
        if (value == "brute_force") {
          measurementSettings.dipoleInteractionHandling =
              DipoleInteractions::kBruteForce;
        } else if (value == "macrocell") {
          measurementSettings.dipoleInteractionHandling =
              DipoleInteractions::kMacrocellMethod;
        } else {
          measurementSettings.dipoleInteractionHandling =
              DipoleInteractions::kNoInteractions;
        }
      } else if (token == "output") {
        measurementSettings.outputPath = value;
      } else if (token == "structure_path") {
        measurementSettings.structurePath = value;
      } else if (token == "num_particles") {
        measurementSettings.numParticles = std::stod(value);
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
      }
    }
  }
  return measurementSettings;
}
