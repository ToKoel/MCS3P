//
//  ParserTest.cpp
//  MCS3P
//
//  Created by Tobias Köhler on 12.04.22.
//

#include "ParserTest.hpp"
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
