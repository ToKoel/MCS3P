//
//  ExchangeConstants.h
//  MCS3P
//
//  Created by Tobias Köhler on 07.02.22.
//

#ifndef HelperStructs_h
#define HelperStructs_h

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "constants.h"
#include <type_traits>

enum class MeasurementType {
    kNone,
    kMvsT,
    kMvsH,
    kSpinStructure,
    kTest
};

enum class DipoleInteractions{ kNoInteractions=0, kBruteForce, kMacrocellMethod };
static std::vector<std::string> DipoleInteractionTypes =
{
    "noDip",
    "brute",
    "macrocell"
};

enum class StructuralPositions{ kOctahedral, kTetrahedral };
static std::vector<std::string> StructurePositionTypes =
{
    "octahedral",
    "tetrahedral"
};

struct ExchangeConstants{
    double FeTT = 0.0;
    double FeOO = 0.0;
    double FeTO = 0.0;
    double FeOO_APB = 0.0;
};

struct LatticeParameters{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
};

struct Angles{
    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;
};

struct LinalgVector{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    
    double selfDot(){
        return x*x + y*y + z*z;
    }
    
    LinalgVector operator-(const LinalgVector& vec2) const{
        return {x-vec2.x, y-vec2.y, z-vec2.z};
    }
    
    LinalgVector operator-() const{
        return {-x, -y, -z};
    }
    
    LinalgVector operator/(double c) const{
        return {x/c, y/c, z/c};
    }
    
    LinalgVector operator+(const LinalgVector& vec2) const{
        return {x+vec2.x, y+vec2.y, z+vec2.z};
    }
    
    LinalgVector operator*(double c) const{
        return {x*c, y*c, z*c};
    }
    
    LinalgVector& operator+=(const LinalgVector& vec2){
        x += vec2.x;
        y += vec2.y;
        z += vec2.z;
        return *this;
    }
    
    LinalgVector& operator+=(double c){
        x += c;
        y += c;
        z += c;
        return *this;
    }
    
    LinalgVector& operator*= (double c){
        x *= c;
        y *= c;
        z *= c;
        return *this;
    }
    
    double dot(const LinalgVector& vec2) const{
        return x*vec2.x + y*vec2.y + z*vec2.z;
    }
    
    LinalgVector cross(const LinalgVector& vec2) const{
        return {y*vec2.z - z*vec2.y,
                z*vec2.x - x*vec2.z,
                x*vec2.y - y*vec2.x};
    }
    
    void normalize(){
        double length = std::sqrt(x*x + y*y + z*z);
        if(length != 0.0){
            x /= length;
            y /= length;
            z /= length;
        }
    }
    
    void rotate(double alpha, double beta, double gamma, double center = 0.0){
        double alphaRad = alpha * PI/180;
        double betaRad = beta * PI/180;
        double gammaRad = gamma * PI/180;
        
        double temp_x, temp_y, temp_z = 0.0;
        
        if(center != 0.0){
            temp_x = x-center;
            temp_y = y-center;
            temp_z = z-center;
        } else {
            temp_x = x;
            temp_y = y;
            temp_z = z;
        }
        
        // x-rotation with rotation matrix
        x = temp_x*1.0 + temp_y*0.0 + temp_z*0.0;
        y = temp_x*0.0 + temp_y*std::cos(alphaRad) + temp_z*-(std::sin(alphaRad));
        z = temp_x*0.0 + temp_y*std::sin(alphaRad) + temp_z*std::cos(alphaRad);
        
        temp_x = x;
        temp_y = y;
        temp_z = z;
        
        // y-rotation
        x = temp_x*std::cos(betaRad) + temp_y*0.0 + temp_z*std::sin(betaRad);
        y = temp_x*0.0 + temp_y*1.0 + temp_z*0.0;
        z = temp_x*-(std::sin(betaRad)) + temp_y*0.0 + temp_z*std::cos(betaRad);
        
        temp_x = x;
        temp_y = y;
        temp_z = z;

        // z-rotation
        x = temp_x*std::cos(gammaRad) + temp_y*-(std::sin(gammaRad)) + temp_z*0.0;
        y = temp_x*std::sin(gammaRad) + temp_y*std::cos(gammaRad) + temp_z*0.0;
        z = temp_x*0.0 + temp_y*0.0 + temp_z*1.0;

        if(center != 0.0){
            x += center;
            y += center;
            z += center;
        }
    }
};

struct MeasurementSettings{
    MeasurementType measurementType = MeasurementType::kNone;
    std::string outputPath = "";
    std::string structurePath = "";
    int numParticles = 0.0;
    int particleSize = 0.0;
    double measurementField = 0.0;
    double temperature = 0.0;
    DipoleInteractions dipoleInteractionHandling = DipoleInteractions::kNoInteractions;
    double monteCarloSteps = 0.0;
    int averagingSteps = 0;
    int numOrientations = 0;
    double coolingField = 0.0;
    double upperTemperatureLimit = 0.0;
    double lowerTemperatureLimit = 0.0;
    double temperatureStepSize = 0.0;
    Angles angles;
    double macrocellSize = 0.0;
    double anisotropyConstant = 0.0;
    bool ZFC = false;
    bool FC = false;
    bool APB = false;
    double upperFieldLimit = 0.0;
    double lowerFieldLimit = 0.0;
    double fieldStepSize = 0.0;
    double startingTemperature = 0.0;
    double coolingSteps = 0.0;
    double sigma = 0.0;
    LatticeParameters latticeParameters;
    ExchangeConstants exchangeConstants;
    double particleCenter = 0.0;
    int totalNumAtoms = 0;
};

struct StructureProperties{
    int numberOfAtoms;
    std::vector<LinalgVector> positionVectors;
    std::vector<StructuralPositions> positionIDs;
    std::vector<bool> isAPB;
};


#endif
