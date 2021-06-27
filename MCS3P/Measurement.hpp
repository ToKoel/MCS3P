#ifndef Measurement_hpp
#define Measurement_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Crystal.hpp"
#include "randNumGenerator.hpp"
#include "omp.h"
#include "ProgressBar.hpp"



template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    if(start < stop){
        for (T value = start; value < stop; value += step)
        values.push_back(value);
    } else{
        for (T value = start; value > stop; value -= step)
            values.push_back(value);
    }
    return values;
}

inline std::string num_to_string(double number, int precision){
    std::string double_string = std::to_string((int)(number))+"d";
    double multiplier = 10.0;
    int A = 0;
    for(int i=0; i<precision; i++){
        A = (int)(number*multiplier - (int)(number*multiplier/10.0)*10);
        double_string += std::to_string(A);
        multiplier *= 10.0;
    }
    return double_string;
}

class Measurement{
public:
    int particleSize;
    bool antiphaseBoundary;
    double surfaceLayerThickness;
    bool dipoleInteractions;
    int steps;
    bool vacancies;
    bool cutOff;
    std::string fileName;
    
    Measurement(int particleSize, bool antiphaseBoundary, double surfaceLayerThickness, bool dipoleInteractions,bool vacancies, bool cutOff, int steps);
};

class MvsTMeasurement{
public:
    // general settings
    std::string dipoleInteractions;
    double steps;
    int averaging_steps;
    int numOrientations;
    double measurement_field;
    double cooling_field;
    double macrocell_size;

    // temperature sweep settings
    double TUpperLimit;
    double TLowerLimit;
    double TstepSize;
    
    int TnumberOfSteps;
    std::vector<double> Tdown;
    std::vector<double> Tup;
    
    std::vector<double> mxZFC;
    std::vector<double> myZFC;
    std::vector<double> mzZFC;
    
    std::vector<double> mxCFZ;
    std::vector<double> myCFZ;
    std::vector<double> mzCFZ;
    
    std::vector<double> mxFC;
    std::vector<double> myFC;
    std::vector<double> mzFC;
    
    // class initializer
    MvsTMeasurement(std::string dipoleInteractions, double steps,
                    int averaging_steps, int numOrientations,
                    double measurement_field, double cooling_field,
                    double TUpperLimit, double TLowerLimit,
                    double TstepSize, double macrocell_size);
    
    // function to run the simulation
    void temperatureSweep(std::string output_dir, std::string structure_filename,
                          double FeTT, double FeOO, double FeTO, double FeOO_APB,
                          double anisotropyConstant,bool ZFC, bool FC, double center,
                          double lattice_a, double lattice_b, double lattice_c, double sigma);
};


// ------- field sweep class ----------------------------
class MvsBMeasurement{
public:
    // general settings
    std::string dipoleInteractions;
    int steps;
    int numOrientations;
    double temperature;
    double macrocell_size;
    
    // cooling setup
    double startTemp;
    double tempStep;
    double coolingField;
    int coolingSteps;
    
    // field sweep settings
    int BnumberOfSteps;
    double BUpperLimit;
    double BLowerLimit;
    double BstepSize;
    std::vector<double> magneticField;
    
    // MvsB arrays for the x, y and z components
    std::vector<double> mX;
    std::vector<double> mY;
    std::vector<double> mZ;
    
    // class initializer
    MvsBMeasurement(std::string dipoleInteractions,
                    int steps, int numOrientations,
                    double temperature, double BUpperLimit,
                    double BLowerLimit, double BstepSize,
                    double startTemp, double tempStep,
                    double coolingField, int coolingSteps, double macrocell_size);
    
    // function to run the simulation
    void fieldSweep(std::string output_dir, std::string structure_filename,
                    double FeTT, double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant,
                    double lattice_a, double lattice_b, double lattice_c,double center, double sigma);
};
// ---------------------------------------------------------
// ------- spin structure class ----------------------------
class spinStructure{
public:
    std::string output_dir;
    std::string structure_filename;
    
    std::string dipoleInteractions;
    int steps;
    double magneticField;
    double temperature;
    double macrocell_size;
    
    spinStructure(std::string dipoleInteractions,int steps, double magneticField, double temperature, double macrocell_size, std::string output_dir, std::string structure_filename);
    
    void spinStructureMeasurement(double FeTT,double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant, double alpha,double beta , double gamma, double center,double lattice_a, double lattice_b, double lattice_c, double sigma);
};
// ---------------------------------------------------------

// wrapper function for M vs. B measurements
void run_MvsB(std::string output_dir, std::string structure_filename,std::string dipoleInteractions, int steps, int numOrientations, double temperature, double BUpperLimit, double BLowerLimit, double BstepSize, double coolingField, int coolingSteps, double startTemp, double tempStep, double FeTT, double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant, double macrocell_size, double center, double lattice_a, double lattice_b, double lattice_c, double sigma);

// wrapper function for M vs. T measurements
void run_MvsT(std::string output_dir, std::string structure_filename,std::string dipoleInteractions,
              double steps,int averaging_steps, int numOrientations, double measurement_field,
              double cooling_field, double TUpperLimit, double TLowerLimit, double TstepSize,
              double FeTT, double FeOO, double FeTO, double FeOO_APB, double anisotropyConstant,
              bool ZFC, bool FC, double macrocell_size, double center,double lattice_a, double lattice_b,
              double lattice_c, double sigma);

// wrapper function for spin structure calculations
void run_spinstructure(std::string dipoleInteractions,int steps, double magneticField, double temperature, std::string output_dir, std::string structure_filename,double FeTT,double FeOO, double FeTO, double FeOO_APB,double anisotropyConstant, double alpha,double beta, double gamma , double macrocell_size, double center, double lattice_a, double lattice_b, double lattice_c, double sigma);



#endif /* Measurement_hpp */
