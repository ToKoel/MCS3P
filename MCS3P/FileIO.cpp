//
//  FileIO.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 23.07.22.
//
#include "FileIO.hpp"

namespace fileIO {

FileHandler::FileHandler(utility::MeasurementSettings measurementSettings): measurementSettings(measurementSettings){};

std::string FileHandler::generateOutputHeader() {
  std::stringstream stream;
  switch (measurementSettings.measurementType) {
      case utility::MeasurementType::kMvsT:
    stream << "# structure file: " << measurementSettings.structureFile << "\n";
    stream << "# dipole interactions: "
           << utility::DipoleInteractionTypes[static_cast<size_t>(
                  measurementSettings.dipoleInteractionHandling)]
           << "\n";
    stream << "# steps: " << measurementSettings.monteCarloSteps << "\n";
    stream << "# number of orientations: "
           << measurementSettings.numOrientations << "\n";
    stream << "# measurement field: " << measurementSettings.measurementField
           << "\n";
    stream << "# cooling field: " << measurementSettings.coolingField << "\n";
    stream << "# upper temperature limit: "
           << measurementSettings.upperTemperatureLimit << "\n";
    stream << "# lower temperature limit: "
           << measurementSettings.lowerTemperatureLimit << "\n";
    stream << "# temperature step size: "
           << measurementSettings.temperatureStepSize << "\n";
    stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
    stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
    stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
    stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB
           << "\n";
    stream << "# lattice parameters: "
           << measurementSettings.latticeParameters.a;
    stream << ", " << measurementSettings.latticeParameters.b;
    stream << ", " << measurementSettings.latticeParameters.c << "\n";
    stream << "# sigma: " << measurementSettings.sigma << "\n";
    stream << "# "
           << "\n";
    stream << "# field (T)      M_x(B)      M_y(B)      M_z(B)"
           << "\n";
    break;

      case utility::MeasurementType::kMvsH:
          stream << "# structure file: " << measurementSettings.structureFile << "\n";
    stream << "# dipole interactions: "
           << utility::DipoleInteractionTypes[static_cast<size_t>(
                  measurementSettings.dipoleInteractionHandling)]
           << "\n";
    stream << "# steps: " << measurementSettings.monteCarloSteps << "\n";
    stream << "# number of orientations: "
           << measurementSettings.numOrientations << "\n";
    stream << "# temperature: " << measurementSettings.temperature << "\n";
    stream << "# B_upper: " << measurementSettings.upperFieldLimit << "\n";
    stream << "# B_lower: " << measurementSettings.lowerFieldLimit << "\n";
    stream << "# B_step: " << measurementSettings.fieldStepSize << "\n";
    stream << "# cooling_field: " << measurementSettings.coolingField << "\n";
    stream << "# start_temperature: " << measurementSettings.startingTemperature
           << "\n";
    stream << "# temperature_step: " << measurementSettings.temperatureStepSize
           << "\n";
    stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
    stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
    stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
    stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB
           << "\n";
    stream << "# lattice parameters: "
           << measurementSettings.latticeParameters.a;
    stream << ", " << measurementSettings.latticeParameters.b;
    stream << ", " << measurementSettings.latticeParameters.c << "\n";
    // stream << "# total number of iron atoms: " << totalNumAtoms << std::endl;
    stream << "# " << std::endl;
    stream << "# field (T)      M_x(B)      M_y(B)      M_z(B)"
           << "\n";
    break;

      case utility::MeasurementType::kSpinStructure:
          stream << "# structure_file: " << measurementSettings.structureFile << "\n";
    stream << "# dipole_interactions: "
           << utility::DipoleInteractionTypes[static_cast<size_t>(
                  measurementSettings.dipoleInteractionHandling)]
           << "\n";
    if (measurementSettings.dipoleInteractionHandling ==
        utility::DipoleInteractions::kMacrocellMethod) {
      stream << "# macrocell_size: " << measurementSettings.macrocellSize
             << "\n";
    }
    stream << "# steps: " << measurementSettings.monteCarloSteps << std::endl;
    stream << "# Number of atoms: " << measurementSettings.totalNumAtoms
           << "\n";
    stream << "# temperature: " << measurementSettings.temperature << "\n";
    stream << "# field: " << measurementSettings.measurementField << "\n";
    stream << "# FeTT: " << measurementSettings.exchangeConstants.FeTT << "\n";
    stream << "# FeTO: " << measurementSettings.exchangeConstants.FeTO << "\n";
    stream << "# FeOO: " << measurementSettings.exchangeConstants.FeOO << "\n";
    stream << "# FeOO_APB: " << measurementSettings.exchangeConstants.FeOO_APB
           << "\n";
    stream << "# anisotropy constant: "
           << measurementSettings.anisotropyConstant << "\n";
    stream << "# lattice parameters: "
           << measurementSettings.latticeParameters.a;
    stream << ", " << measurementSettings.latticeParameters.b;
    stream << ", " << measurementSettings.latticeParameters.c << "\n";
    stream << "# Orientation: " << measurementSettings.angles.alpha << ", ";
    stream << measurementSettings.angles.beta << ", "
           << measurementSettings.angles.gamma << "\n";
    stream << "# "
           << "\n";
    stream << "# x      y        z      spinx    spiny    spinz  pos  APB"
           << "\n";
    break;

      case utility::MeasurementType::kNone:
    break;

      case utility::MeasurementType::kTest:
    break;
  }

  return stream.str();
}

std::string FileHandler::generateOutputFilename() {
    std::string structure_path = measurementSettings.structureFile;
    std::string outputFilename = measurementSettings.outputPath + "/" + structure_path.substr(structure_path.find_last_of("/") + 1);
  switch (measurementSettings.measurementType) {
      case utility::MeasurementType::kMvsT:
    outputFilename +=
        "_MvsT_sim_cF" + utility::to_string(measurementSettings.coolingField, 2);
    outputFilename +=
        "T_mF" + utility::to_string(measurementSettings.measurementField, 3);
    outputFilename +=
        "T_" + utility::to_string(measurementSettings.monteCarloSteps, 1) + "steps";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.averagingSteps, 1) + "avsteps";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.numOrientations, 1) + "or";
    outputFilename += "_" + utility::to_string(measurementSettings.sigma, 2) + "sig";
    break;

      case utility::MeasurementType::kMvsH:
    outputFilename +=
        "_Hysteresis_sim_" + utility::to_string(measurementSettings.temperature, 1);
    outputFilename +=
        "K_" + utility::to_string(measurementSettings.monteCarloSteps, 1) + "steps";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.numOrientations, 1) + "or";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.coolingField, 2) + "T";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.coolingSteps, 2) + "cs";
    outputFilename +=
        "_sT" + utility::to_string(measurementSettings.startingTemperature, 1) + "K";
    outputFilename +=
        "_dip" + utility::DipoleInteractionTypes[static_cast<size_t>(
                     measurementSettings.dipoleInteractionHandling)];
    break;

      case utility::MeasurementType::kSpinStructure:
    outputFilename +=
        "_spin_structure_" + utility::to_string(measurementSettings.temperature, 1);
    outputFilename +=
        "K_" + utility::to_string(measurementSettings.monteCarloSteps, 1) + "MCS";
    outputFilename +=
        "_" + utility::to_string(measurementSettings.measurementField, 1) + "T";
    outputFilename +=
        "_dip" + utility::DipoleInteractionTypes[static_cast<size_t>(
                     measurementSettings.dipoleInteractionHandling)];
    break;

      case utility::MeasurementType::kNone:
    break;

      case utility::MeasurementType::kTest:
    break;
  }

  if (measurementSettings.dipoleInteractionHandling ==
      utility::DipoleInteractions::kMacrocellMethod) {
    outputFilename +=
        "_" + utility::to_string(measurementSettings.macrocellSize, 4) + "mcsize";
  }
  return outputFilename;
}

void FileHandler::writeSweepResultsToFile(
                             utility::MeasurementMode mode, std::vector<double> &variableArray,
    std::vector<utility::LinalgVector> &magnetizationArray) {
  std::ofstream sweep;
  switch (mode) {
      case utility::MeasurementMode::ZFC:
    sweep.open(generateOutputFilename() + "_ZFC.txt", std::fstream::out);
    break;
      case utility::MeasurementMode::FC:
    sweep.open(generateOutputFilename() + "_FC.txt", std::fstream::out);
    break;
      case utility::MeasurementMode::Field:
    sweep.open(generateOutputFilename() + ".txt", std::fstream::out);
    break;
  }
  sweep << generateOutputHeader();
  for (size_t currentStep = 0; currentStep < variableArray.size(); currentStep++) {
    sweep << variableArray[currentStep] << " "
          << magnetizationArray[currentStep].x << " "
          << magnetizationArray[currentStep].y << " "
          << magnetizationArray[currentStep].z << "\n";
  }
  sweep.close();
}

void FileHandler::writeSpinStructureResultsToFile(Crystal &crystal) {
  std::ofstream structure;
  structure.open(generateOutputFilename() + ".txt", std::fstream::out);
  structure << generateOutputHeader();
  for (auto atom : crystal.atoms) {
    structure
        << atom.positionVector.x << ", " << atom.positionVector.y << ", "
        << atom.positionVector.z << ", " << atom.spinVector.x << ", "
        << atom.spinVector.y << ", " << atom.spinVector.z << ", "
        << utility::StructurePositionTypes[static_cast<size_t>(atom.structuralPositionID)]
        << ", " << atom.isApbAtom << "\n";
  }
  structure.close();
}

}
