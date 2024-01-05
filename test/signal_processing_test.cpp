#include "hampel.hpp"
#include <gtest/gtest.h>
#include <rbm/rbm.hpp>

using namespace gpw::sp::filter;
using namespace gpw::geometry;

TEST (Hampel, Median) {
    std::vector<int> x {2, 4, 6, 8, 1, 3, 5, 7, 9};
    EXPECT_EQ (median (x.begin(), x.end()), 5);

    std::vector<vec3> v1 {
        vec3 {1, 2, 3},
         vec3 {2, 3, 4},
         vec3 {3, 4, 5}
    };
    auto mid = hampel_filter (v1.begin(), v1.end());
    EXPECT_EQ (mid (1), 2);
    EXPECT_EQ (mid (2), 3);
    EXPECT_EQ (mid (3), 4);

    std::vector<vec3> v2 {
        vec3 { 1, 2, 3}, // norm = 3.7417
        vec3 { 2, 2, 1}, // norm = 3
        vec3 {52, 2, 2}, // norm = 52.0769
        vec3 { 3, 1, 5}, // norm = 5.9161
        vec3 { 2, 1, 2}  // norm = 3
    };
    auto mid1 = hampel_filter (v2.begin(), v2.end());
    EXPECT_EQ (mid1 (1), 1);
    EXPECT_EQ (mid1 (2), 2);
    EXPECT_EQ (mid1 (3), 3);
}

TEST (Hampel, Window3_Case_001) {
    // Arrange.
    std::vector<int> x {0, 1, 2};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end()), 1);
}

// The following four tests based on the following material:
// https://www.mathworks.com/help/dsp/ref/hampelfilter.html .
TEST (Hampel, Window5_Case_002) {
    // Arrange.
    std::vector<int> x {1, 2, 4, 9, 23};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end()), 4);
}

TEST (Hampel, Window5_Case_003) {
    // Arrange.
    std::vector<int> x {2, 4, 9, 23, 8};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end()), 9);
}

TEST (hampel, Window5_Case_004) {
    // Arrange.
    std::vector<int> x {4, 9, 23, 8, 12};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end(), 5), 23);
}

TEST (Hampel, Window5_Case_005) {
    // Arrange.
    std::vector<int> x {4, 9, 23, 8, 12};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end(), 2), 9);
}

TEST (Hampel, Window3_Case_006) {
    // Arrange.
    std::vector<int> x {0, 1, 2};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end(), 0), 1);
}

TEST (Hampel, Window3_Case_007) {
    // Arrange.
    std::vector<int> x {3, 1, 2};

    // Act. Assert.
    EXPECT_EQ (hampel_filter (x.begin(), x.end(), 0), 2);
    EXPECT_EQ (hampel_filter (x.begin(), x.end(), 1), 1);
}

TEST (Hampel, Vec3) {
    std::vector<vec3> v {
        vec3 {1, 2, 3},
         vec3 {2, 3, 4},
         vec3 {3, 4, 5}
    };
    auto mid = hampel_filter (v.begin(), v.end());
    EXPECT_EQ (mid (1), 2);
    EXPECT_EQ (mid (2), 3);
    EXPECT_EQ (mid (3), 4);

    std::vector<vec3> v1 {
        vec3 { 1, 2, 3}, // norm = 3.7417
        vec3 { 2, 2, 1}, // norm = 3
        vec3 {52, 2, 2}, // norm = 52.0769
        vec3 { 3, 1, 5}, // norm = 5.9161
        vec3 { 2, 1, 2}  // norm = 3
    };
    auto mid1 = hampel_filter (v1.begin(), v1.end());
    EXPECT_EQ (mid1 (1), 1);
    EXPECT_EQ (mid1 (2), 2);
    EXPECT_EQ (mid1 (3), 3);
}
