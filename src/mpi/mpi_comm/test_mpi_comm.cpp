
#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "mpi_headers/mpi_comm.h"
#include "mpi_headers/mpi_test_env.h"
#include "mpi_headers/mpi_utils.h"

class MPI_send_pos_vel_test : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;

        double *vxbuf;
        double *vybuf;
        double *vzbuf;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                ASSERT_EQ(nprocs, 4);

                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                sys->natoms = 4;

                sys->vx = new double[1]();
                sys->vy = new double[1]();
                sys->vz = new double[1]();

                if (proc_id == 0) {

                        sys->rx = new double[4]{0, 10, 20, 30};
                        sys->ry = new double[4]{0, 10, 20, 30};
                        sys->rz = new double[4]{0, 10, 20, 30};

                        vxbuf = new double[4]{0, 1, 2, 3};
                        vybuf = new double[4]{0, 1, 2, 3};
                        vzbuf = new double[4]{0, 1, 2, 3};
                } else {

                        sys->rx = new double[4];
                        sys->ry = new double[4];
                        sys->rz = new double[4];

                        vxbuf = nullptr;
                        vybuf = nullptr;
                        vzbuf = nullptr;
                }
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;

                delete[] sys->vx;
                delete[] sys->vy;
                delete[] sys->vz;

                if (proc_id == 0) {
                        delete[] vxbuf;
                        delete[] vybuf;
                        delete[] vzbuf;
                }

                delete sys;
        }
};

TEST_F(MPI_send_pos_vel_test, rbip) {

        arr_seg_t proc_seg;

        init_segments(nprocs, proc_id, &proc_seg, sys->natoms);

        mpi_send_pos_vel(nprocs, &proc_seg, sys, vxbuf, vybuf, vzbuf);

        EXPECT_DOUBLE_EQ(sys->rx[proc_id], 10 * proc_id);
        EXPECT_DOUBLE_EQ(sys->ry[proc_id], 10 * proc_id);
        EXPECT_DOUBLE_EQ(sys->rz[proc_id], 10 * proc_id);

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

int main(int argc, char **argv) {

        std::string command_line_arg(argc == 2 ? argv[1] : "");
        testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new MPITestEnv(command_line_arg));

        return RUN_ALL_TESTS();
}