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

class test_init_segments : public ::testing::Test {

      protected:
        int size, nprocs;
        arr_seg_t proc_seg;

        virtual void SetUp() override {
                size = 998;
                nprocs = 4;
        }

        virtual void TearDown() override {}
};

TEST_F(test_init_segments, all) {

        int idx_arr[4]{0, 250, 500, 749};
        int size_arr[4]{250, 250, 249, 249};

        for (int proc_id = 0; proc_id < nprocs; ++proc_id) {
                init_segments(nprocs, proc_id, &proc_seg, size);

                EXPECT_EQ(proc_seg.id, proc_id);
                EXPECT_EQ(proc_seg.idx, idx_arr[proc_id]);
                EXPECT_EQ(proc_seg.size, size_arr[proc_id]);
        }
}
