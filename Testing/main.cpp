//main.cpp
namespace testing {
namespace internal {
using UInt64 = uint64_t;
using Int64 = int64_t;
using Int32 = int32_t;
} // namespace internal
} // namespace testing

#include "gtest/gtest.h"

class Adder{
public:
    int add(int a, int b){
        return a+b;
    }
};

class AddingTest : public testing::Test{
public:
    Adder adder;
}


TEST_F(AddingTest, test1){
    ASSERT_EQ(2, adder.add(1,1));
}

int main(int argc, const char * argv[]) {
    ::testing::InitGoogleTest(&argc, (char**)argv);
    return RUN_ALL_TESTS();
}
