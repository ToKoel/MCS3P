#ifndef Crystal_h
#define Crystal_h

#include <stdio.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include "Atom.hpp"
#include "randNumGenerator.hpp"
#include "constants.hpp"
#include "HelperStructs.hpp"
#include "Parsers.hpp"

class Crystal{
public:
    std::vector<Atom> atoms;
    utility::MeasurementSettings measurementSettings;
    Crystal(utility::MeasurementSettings);
    Crystal(){};
    
    void rotateCrystal();
    void randomOrientation();
    void alignAlongRandomVector();
    void structureSnapshot(std::string filename);
    void resetStructure();
    int outputStats();
    void setSigma();
    void generateDipoleLists();
    void generateNeighbourLists();
    void initializeStructureFromFile();
    int getNumberOfAtoms();
    //std::vector<Macrocell> macrocells;
    void generateMacrocells();
    void saveMacrocells(std::string filename);
    void performMonteCarloSteps(double numberOfSteps, utility::Environment environment);
    utility::LinalgVector getNetSpinVector();
};

#endif
