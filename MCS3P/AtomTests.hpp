//
//  AtomTest.h
//  MCS3P
//
//  Created by Tobias Köhler on 13.04.22.
//

#ifndef AtomTest_h
#define AtomTest_h

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Atom.hpp"

class AtomTest : public ::testing::Test {
protected:
    Atom* atom;
    Atom* atom2;
    Atom* atom3;
    Atom* atom4;
    
    void SetUp(){
        MeasurementSettings measurementSettings;
        measurementSettings.exchangeConstants = {-21.00, -8.60, -28.10, -106.28};
        measurementSettings.anisotropyConstant = 3.25e-25;
        measurementSettings.latticeParameters.a = 8.395;
        measurementSettings.latticeParameters.b = 8.395;
        measurementSettings.latticeParameters.c = 8.395;
        
        atom = new Atom(measurementSettings,
                        {0.2, 0.3, 0.4},
                        StructuralPositions::kOctahedral,
                        true);
        atom->spinVector = {0.42426407, 0.56568542, 0.70710678};
        
        atom2 = new Atom(measurementSettings,
                         {0.0, 0.0, 0.0},
                         StructuralPositions::kOctahedral,
                         true);
        atom2->spinVector = {0.0       , 0.9701425 , 0.24253563};
        
        atom3 = new Atom(measurementSettings,
                         {0.1, 0.1, 0.2},
                         StructuralPositions::kOctahedral,
                         false);
        atom3->spinVector = {0.13018891, 0.39056673, 0.91132238};
        
        atom4 = new Atom(measurementSettings,
                         {0.2, 0.1, 0.2},
                         StructuralPositions::kTetrahedral,
                         false);
        atom4->spinVector = {0.87287156, 0.21821789, 0.43643578};
        
        atom->neighboursParallel.push_back(atom2);
        atom->neighboursParallel.push_back(atom3);
        atom->neighboursAntiparallel.push_back(atom4);
        
        atom->allOtherAtomsInCrystal.push_back(atom2);
        atom->allOtherAtomsInCrystal.push_back(atom3);
        atom->allOtherAtomsInCrystal.push_back(atom4);
        
        atom->inv_distances_cubed.push_back(0.010822831617974598);
        atom->inv_distances_cubed.push_back(0.06259997134581934);
        atom->inv_distances_cubed.push_back(0.0746969584062022);
        atom->distanceVectors.push_back({-0.3713906763541038, -0.5570860145311556, -0.7427813527082076});
        atom->distanceVectors.push_back({-0.33333333333333337, -0.6666666666666666, -0.6666666666666667});
        atom->distanceVectors.push_back({0.0        , -0.7071067811865475, -0.7071067811865476});
    }
    
    void TearDown(){
        delete atom;
        delete atom2;
        delete atom3;
        delete atom4;
    }
};

TEST_F(AtomTest, atomIsInitialized){
    ASSERT_EQ(atom->exchangeConstants.FeTT, -21.0);
    ASSERT_EQ(atom->exchangeConstants.FeTO, -28.10);
    ASSERT_EQ(atom->exchangeConstants.FeOO, -8.60);
    ASSERT_EQ(atom->exchangeConstants.FeOO_APB, -106.28);
    ASSERT_THAT(atom->structuralPositionID, StructuralPositions::kOctahedral);
    ASSERT_EQ(atom->positionVector.x, 0.2);
    ASSERT_EQ(atom->positionVector.y, 0.3);
    ASSERT_EQ(atom->positionVector.z, 0.4);
    ASSERT_TRUE(atom->isApbAtom);
}

TEST_F(AtomTest, atomAnisotropyEnergyIsCalculatedCorrectly){
    ASSERT_DOUBLE_EQ(atom->anisotropy(), -6.252999920891827e-26);
}

TEST_F(AtomTest, atomZeemanEnergyIsCalculatedCorrectly){
    ASSERT_DOUBLE_EQ(atom->zeeman(1.2), -2.360777556624346e-23);
}

TEST_F(AtomTest, atomExchangeEnergyIsCalculatedCorrectly){
    ASSERT_DOUBLE_EQ(atom->exchange(), 1.477523609648334e-21);
}

TEST_F(AtomTest, atomDipoleEnergyIsCalculatedCorrectly){
    ASSERT_DOUBLE_EQ(atom->dipole(), -3.456985088289361e-23);
}

TEST_F(AtomTest, atomSpinIsFlipped){
    atom->spin_flip();
    ASSERT_EQ(atom->spinVector.x, -0.42426407);
}

#endif /* AtomTest_h */
