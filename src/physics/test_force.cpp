#include "gtest/gtest.h"

/* inclusion of source file is needed
 * since we want to test a static function  */
extern "C" {
#include "force.c"
}

TEST(pbc_test, inside) {
        ASSERT_DOUBLE_EQ(pbc(3.0, 10.0 * 0.5), 3.0);
        ASSERT_DOUBLE_EQ(pbc(-2.5, 10.0 * 0.5), -2.5);
}

TEST(pbc_test, outside) {
        ASSERT_DOUBLE_EQ(pbc(7.0, 10.0 * 0.5), -3.0);
        ASSERT_DOUBLE_EQ(pbc(-8.0, 10.0 * 0.5), 2.0);
}

TEST(ForceTestSingle, single) {
        mdsys_t *sys = new mdsys_t;
        sys->natoms = 1;

#if defined(MPI_ENABLED)
        arr_seg_t proc_seg;
        sys->proc_seg = &proc_seg;
        proc_seg.size = sys->natoms;
        proc_seg.idx = 0;
#endif

        sys->rx = new double[1];
        sys->ry = new double[1];
        sys->rz = new double[1];

        sys->fx = new double[1];
        sys->fy = new double[1];
        sys->fz = new double[1];

        force(sys);

        ASSERT_EQ(sys->fx[0], 0.0);
        ASSERT_EQ(sys->fy[0], 0.0);
        ASSERT_EQ(sys->fz[0], 0.0);
        ASSERT_DOUBLE_EQ(sys->epot, 0.0);

        delete[] sys->fx;
        delete[] sys->fy;
        delete[] sys->fz;

#if defined(MPI_ENABLED)
        sys->proc_seg = nullptr;
#endif
        delete sys;
}

class ForceTest : public ::testing::TestWithParam<double> {

      protected:
        mdsys_t *sys;
        int eps_param = GetParam();

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->epsilon = eps_param;
                sys->sigma = 1.0;
                sys->box = 10.0;

#if defined(MPI_ENABLED)
                arr_seg_t proc_seg;
                sys->proc_seg = &proc_seg;
#endif

                sys->rx = new double[2]();
                sys->ry = new double[2]();
                sys->rz = new double[2]();

                sys->fx = new double[2];
                sys->fy = new double[2];
                sys->fz = new double[2];

                sys->rx[0] = -1.0;
                sys->rx[1] = 1.0;

                // forces and will be zeroed by azzero
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;

                delete[] sys->fx;
                delete[] sys->fy;
                delete[] sys->fz;

#if defined(MPI_ENABLED)
                sys->proc_seg = nullptr;
#endif
                delete sys;
        }
};

TEST_P(ForceTest, shortrange) {

        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys->natoms, 2);
        ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
        ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
        ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
        ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
        ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
        ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

#if defined(MPI_ENABLED)
        sys->proc_seg->size = 2;
        sys->proc_seg->idx = 0;
#endif
        sys->rcut = 0.5;

        force(sys);

        EXPECT_DOUBLE_EQ(sys->fx[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->fx[1], 0.0);
        EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
        EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
        EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

TEST_P(ForceTest, longrange) {

        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys->natoms, 2);
        ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
        ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
        ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
        ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
        ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
        ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

#if defined(MPI_ENABLED)
        sys->proc_seg->size = 2;
        sys->proc_seg->idx = 0;
#endif

        sys->rcut = 4.0;

        force(sys);

        double exp_epot = -eps_param * 63.0 / 1024.0;
        double exp_ff = -eps_param * 93.0 / 512.0;

        EXPECT_DOUBLE_EQ(sys->fx[0], -exp_ff);
        EXPECT_DOUBLE_EQ(sys->fx[1], exp_ff);
        EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
        EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
        EXPECT_DOUBLE_EQ(sys->epot, exp_epot);
}

INSTANTIATE_TEST_SUITE_P(ForceTest_parametric, ForceTest,
                         ::testing::Values(0.0, 0.5, 1.0));
