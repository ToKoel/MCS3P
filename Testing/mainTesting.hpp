//
//  mainTesting.h
//  Testing
//
//  Created by Tobias Köhler on 10.04.22.
//

#ifndef mainTesting_hpp
#define mainTesting_hpp

#include "gtest/gtest.h"

int mainTesting(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif /* mainTesting_hpp */
