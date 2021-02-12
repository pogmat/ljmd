
#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "mpi_headers/mpi_comm.h"
#include "mpi_headers/mpi_test_env.h"
#include "mpi_headers/mpi_utils.h"
/*
class MPI_send_pos_vel_test : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;
                sys->nprocs = nprocs;
                sys->proc_id = proc_id;

                arr_seg_t proc_seg;
                sys->proc_seg = &proc_seg;

                if (proc_id == 0) {

                        sys->rx = new double[4]{0, 10, 20, 30};
                        sys->ry = new double[4]{0, 10, 20, 30};
                        sys->rz = new double[4]{0, 10, 20, 30};

                        sys->vx = new double[4]{0, 1, 2, 3};
                        sys->vy = new double[4]{0, 1, 2, 3};
                        sys->vz = new double[4]{0, 1, 2, 3};
                } else {

                        sys->rx = new double[4];
                        sys->ry = new double[4];
                        sys->rz = new double[4];

                        sys->vx = new double[1]();
                        sys->vy = new double[1]();
                        sys->vz = new double[1]();
                }
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;

                delete[] sys->vx;
                delete[] sys->vy;
                delete[] sys->vz;

                delete sys;
        }
};

TEST_F(MPI_send_pos_vel_test, rbip) {

        init_segments(nprocs, proc_id, sys->proc_seg, sys->natoms);

        mpi_send_pos_vel(sys);

        for (int p = 0; p < nprocs; ++p) {
                EXPECT_DOUBLE_EQ(sys->rx[p], 10 * p);
                EXPECT_DOUBLE_EQ(sys->ry[p], 10 * p);
                EXPECT_DOUBLE_EQ(sys->rz[p], 10 * p);
        }

        EXPECT_DOUBLE_EQ(sys->vx[0], proc_id);
        EXPECT_DOUBLE_EQ(sys->vy[0], proc_id);
        EXPECT_DOUBLE_EQ(sys->vz[0], proc_id);
}

class MPI_xch_pos : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;
        arr_seg_t proc_seg;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;

                sys->rx = new double[4];
                sys->ry = new double[4];
                sys->rz = new double[4];

                sys->rx[proc_id] = 10.0 * proc_id;
                sys->ry[proc_id] = 10.0 * proc_id;
                sys->rz[proc_id] = 10.0 * proc_id;
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;

                delete sys;
        }
};

TEST_F(MPI_xch_pos, basic) {

        init_segments(nprocs, proc_id, &proc_seg, sys->natoms);

        int count[nprocs];
        int offsets[nprocs];
        mpi_collective_comm_arrays(nprocs, proc_seg.splitting, count, offsets);

        mpi_exchange_positions(sys, count, offsets);

        for (int p = 0; p < nprocs; ++p) {
                EXPECT_DOUBLE_EQ(sys->rx[p], 10.0 * p);
                EXPECT_DOUBLE_EQ(sys->ry[p], 10.0 * p);
                EXPECT_DOUBLE_EQ(sys->rz[p], 10.0 * p);
        }
}

class MPI_red_UKT : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;

                sys->ekin = 1.0;
                sys->temp = 10.0;
                sys->epot = 2.0;
        }

        void TearDown() { delete sys; }
};

TEST_F(MPI_red_UKT, basic) {

        mpi_reduce_UKT(sys);

        if (proc_id == 0) {
                EXPECT_DOUBLE_EQ(sys->ekin, nprocs);
                EXPECT_DOUBLE_EQ(sys->temp, 10.0 * nprocs);
                EXPECT_DOUBLE_EQ(sys->epot, 2.0 * nprocs);
        }
}

*/

class MPI_broadc_pos_test : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;
                sys->nprocs = nprocs;
                sys->proc_id = proc_id;

                if (proc_id == 0) {

                        sys->rx = new double[4]{0, 10, 20, 30};
                        sys->ry = new double[4]{0, 10, 20, 30};
                        sys->rz = new double[4]{0, 10, 20, 30};

                } else {

                        sys->rx = new double[4];
                        sys->ry = new double[4];
                        sys->rz = new double[4];
                }
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;

                delete sys;
        }
};

TEST_F(MPI_broadc_pos_test, simple) {

        mpi_broadcast_pos(sys);

        for (int p = 0; p < sys->natoms; ++p) {
                EXPECT_DOUBLE_EQ(sys->rx[p], 10 * p);
                EXPECT_DOUBLE_EQ(sys->ry[p], 10 * p);
                EXPECT_DOUBLE_EQ(sys->rz[p], 10 * p);
        }
}

class MPI_reduce_forces_u_test : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;
                sys->nprocs = nprocs;
                sys->proc_id = proc_id;

                sys->fx = new double[4]();
                sys->fy = new double[4]();
                sys->fz = new double[4]();

                sys->fx[proc_id] = proc_id;
                sys->fy[proc_id] = proc_id;
                sys->fz[proc_id] = proc_id;

                sys->epot = 10.0;
        }

        void TearDown() {
                delete[] sys->fx;
                delete[] sys->fy;
                delete[] sys->fz;

                delete sys;
        }
};

TEST_F(MPI_reduce_forces_u_test, basic) {

        mpi_reduce_forces(sys);

        if (proc_id == 0) {
                for (int p = 0; p < sys->natoms; ++p) {
                        EXPECT_DOUBLE_EQ(sys->fx[p], p);
                        EXPECT_DOUBLE_EQ(sys->fy[p], p);
                        EXPECT_DOUBLE_EQ(sys->fz[p], p);
                }
        }
        EXPECT_DOUBLE_EQ(sys->epot, 10.0 * nprocs);
}

int main(int argc, char **argv) {

        std::string command_line_arg(argc == 2 ? argv[1] : "");
        testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new MPITestEnv(command_line_arg));

        return RUN_ALL_TESTS();
}