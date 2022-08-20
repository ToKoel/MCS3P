//
//  mainTesting.cpp
//  Testing
//
//  Created by Tobias Köhler on 10.04.22.
//

#include "gtest/gtest.h"

class Adder{
public:
    int add(int a, int b){
        return a+b;
    }
};

class AdderTest : public testing::Test{
public:
    Adder adder;
};

TEST_F(AdderTest, addTwoNumbers){
    ASSERT_EQ(2, adder.add(1,1));
}

//run all tests
int mainTesting(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
