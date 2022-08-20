#include "Measurement.hpp"

void Measurement::temperatureSweep(
    Crystal &crystal, std::vector<double> &temperatureArray,
                                   std::vector<utility::LinalgVector> &magnetizationArray) {
  auto totalNumAtoms = crystal.atoms.size();
  for (size_t i = 0; i < temperatureArray.size(); i++) {
    crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps,
                                   {measurementSettings.measurementField,
        temperatureArray[i]});

    crystal.performMonteCarloSteps(measurementSettings.averagingSteps,
                                   {measurementSettings.measurementField,
        temperatureArray[i]});

    magnetizationArray[i] +=
        crystal.getNetSpinVector() /
        (static_cast<int>(totalNumAtoms) * measurementSettings.numOrientations);
  }
}

void Measurement::fieldSweep(Crystal &crystal,
                             std::vector<double> &temperatureArray,
                             std::vector<double> &fieldArray,
                             std::vector<utility::LinalgVector> &magnetizationArray) {
  auto totalNumAtoms = crystal.atoms.size();

  for (size_t i = 0; i < temperatureArray.size(); i++) {
    crystal.performMonteCarloSteps(
        measurementSettings.coolingSteps * static_cast<double>(totalNumAtoms),
                                   {measurementSettings.coolingField, temperatureArray[i]});
  }

  for (size_t i = 0; i < fieldArray.size(); i++) {
    crystal.performMonteCarloSteps(
                                   measurementSettings.coolingSteps * static_cast<double>(totalNumAtoms), {fieldArray[i],
                                       measurementSettings.temperature});

    magnetizationArray[i] +=
        crystal.getNetSpinVector() /
        (static_cast<int>(totalNumAtoms) *
                             measurementSettings.numOrientations);
  }
}

void Measurement::run_MvsT() {
  std::cout << "\nOutput filename: " << fileHandler.generateOutputFilename() << "\n";
  std::cout << "\n ------- starting M vs T measurement -------\n"
            << "\n";

  auto temperatureStepsCooling =
      utility::arange<double>(measurementSettings.upperTemperatureLimit,
                     measurementSettings.lowerTemperatureLimit,
                     measurementSettings.temperatureStepSize);
  auto temperatureStepsHeating =
      utility::arange<double>(measurementSettings.lowerTemperatureLimit,
                     measurementSettings.upperTemperatureLimit,
                     measurementSettings.temperatureStepSize);

    std::vector<utility::LinalgVector> magnetizationZFC;
    std::vector<utility::LinalgVector> magnetizationFC;

  magnetizationZFC.resize(temperatureStepsHeating.size());
  magnetizationFC.resize(temperatureStepsCooling.size());

  ProgressBar bar;

  Crystal crystal(measurementSettings);

  for (int i = 0; i < measurementSettings.numOrientations; i++) {
    bar.step(static_cast<double>(i) / measurementSettings.numOrientations,
             "orientation " + std::to_string(i + 1));

    crystal.randomOrientation();
    crystal.alignAlongRandomVector();

    if (measurementSettings.ZFC) {
      temperatureSweep(crystal, temperatureStepsHeating, magnetizationZFC);
        fileHandler.writeSweepResultsToFile(utility::MeasurementMode::ZFC, temperatureStepsHeating,
                              magnetizationZFC);
    }

    if (measurementSettings.FC) {
      temperatureSweep(crystal, temperatureStepsCooling, magnetizationFC);
        fileHandler.writeSweepResultsToFile(utility::MeasurementMode::FC, temperatureStepsCooling,
                              magnetizationFC);
    }
  }
}

void Measurement::run_MvsB() {
  std::string outputFileName = fileHandler.generateOutputFilename();
  std::cout << "\nOutput filename: " << outputFileName << "\n";
  std::cout << "\n ------- starting M vs B measurement -------\n"
            << "\n";

  std::vector<double> fieldArray =
    utility::arange<double>(0, measurementSettings.upperFieldLimit,
                     measurementSettings.fieldStepSize);
    std::vector<double> fieldArrayUp = utility::arange<double>(
      measurementSettings.lowerFieldLimit, measurementSettings.upperFieldLimit,
      measurementSettings.fieldStepSize);
    std::vector<double> fieldArrayDown = utility::arange<double>(
      measurementSettings.upperFieldLimit, measurementSettings.lowerFieldLimit,
      measurementSettings.fieldStepSize);
    std::vector<double> temperatureArray = utility::arange<double>(
      measurementSettings.startingTemperature, measurementSettings.temperature,
      measurementSettings.temperatureStepSize);
  fieldArray.insert(fieldArray.end(), fieldArrayDown.begin(),
                    fieldArrayDown.end());
  fieldArray.insert(fieldArray.end(), fieldArrayUp.begin(), fieldArrayUp.end());
  fieldArray.push_back(measurementSettings.upperFieldLimit);
    std::vector<utility::LinalgVector> magnetization;
  magnetization.resize(fieldArray.size());

  ProgressBar bar;
  for (int i = 0; i < measurementSettings.numOrientations; i++) {
    bar.step(static_cast<double>(i) / measurementSettings.numOrientations,
             "orientation " + std::to_string(i + 1));
    Crystal crystal(measurementSettings);
    crystal.randomOrientation();
    fieldSweep(crystal, temperatureArray, fieldArray, magnetization);
  }
    fileHandler.writeSweepResultsToFile(utility::MeasurementMode::Field, fieldArray, magnetization);
}

void Measurement::run_spinstructure() {
  std::cout << "\n ------- starting spin structure measurement -------\n"
            << "\n";

  Crystal crystal(measurementSettings);
    std::cout << "Crystal initialized" << std::endl;

  // simulated annealing
    if (!(measurementSettings.startingTemperature < measurementSettings.temperature)){
  std::vector<double> coolingArray =
    utility::arange(measurementSettings.startingTemperature,
             measurementSettings.temperature, measurementSettings.coolingSteps);
  coolingArray.push_back(measurementSettings.temperature);
    std::cout << "starting cooling" << std::endl;
  for (double temperatureStep : coolingArray) {
    crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps / 100,
                                   {measurementSettings.measurementField,
        temperatureStep});
  }
    }
    std::cout << "starting relaxation" << std::endl;
  // relaxation
    utility::Environment environment{measurementSettings.measurementField, measurementSettings.temperature};
  crystal.performMonteCarloSteps(measurementSettings.monteCarloSteps, environment);

    fileHandler.writeSpinStructureResultsToFile(crystal);
}
