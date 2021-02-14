#include "common.h"
#include "gtest/gtest.h"

#include <time.h>

class azzero_test : public ::testing::TestWithParam<int> {

      protected:
        double *buf;
        int size = GetParam();

        virtual void SetUp() override {
                buf = new double[size];
                for (int i = 0; i < size; ++i)
                        buf[i] = static_cast<double>(i + 1);
        }

        virtual void TearDown() override { delete[] buf; }
};

TEST(test_wallclock, one_second) {
        struct timespec t;
        t.tv_sec = 0;
        t.tv_nsec = 100000000;

        double t_test = -wallclock();
        nanosleep(&t, NULL);
        t_test += wallclock();
        ASSERT_TRUE((t_test > 0.09));
        ASSERT_TRUE((t_test < 0.11));
}

TEST_P(azzero_test, fill_empty) {

        for (int i = 0; i < size; ++i)
                ASSERT_DOUBLE_EQ(buf[i], static_cast<double>(i + 1));
        azzero(buf, size);
        for (int i = 0; i < size; ++i)
                ASSERT_DOUBLE_EQ(buf[i], 0.0);
}

INSTANTIATE_TEST_SUITE_P(azzero_test_parametric, azzero_test,
                         ::testing::Values(0, 1, 10, 50));
