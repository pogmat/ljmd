#include "gtest/gtest.h"
#include <math.h>

#include "mpi_headers/mpi_utils.h"

class test_split_dimension : public ::testing::TestWithParam<int> {

      protected:
        int size, nprocs;
        int splitting[3];

        virtual void SetUp() override {
                size = 1000;
                nprocs = GetParam();
                ASSERT_TRUE(nprocs > 0);
        }

        virtual void TearDown() override {}
};

TEST_P(test_split_dimension, all) {
        split_dimension(nprocs, size, splitting);

        float base = ((float)size) / ((float)nprocs);
        int h1 = (int)ceil(base);
        int h2 = (int)floor(base);

        ASSERT_TRUE((splitting[0] > 0));
        EXPECT_TRUE((splitting[0] <= nprocs));

        EXPECT_TRUE((splitting[1] - splitting[2]) >= 0);
        EXPECT_TRUE((splitting[1] - splitting[2]) <= 1);

        EXPECT_EQ(splitting[0] * splitting[1] +
                      (nprocs - splitting[0]) * splitting[2],
                  size);
}

INSTANTIATE_TEST_SUITE_P(test_split_dimension_p, test_split_dimension,
                         ::testing::Values(1, 4, 10));