
#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "mpi_headers/mpi_comm.h"
#include "mpi_headers/mpi_test_env.h"
#include "mpi_headers/mpi_utils.h"

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
				
				sys->r = new vec3_t[4];

                if (proc_id == 0) {

						
					
						for (int i=0; i < 4; ++i ) {
							sys->r[i].x = i*10.0;
							sys->r[i].y = i*10.0;
							sys->r[i].z = i*10.0;
						}
                       

                } 
        }

        void TearDown() {
                delete[] sys->r;

                delete sys;
        }
};

TEST_F(MPI_broadc_pos_test, simple) {

        mpi_broadcast_pos(sys);

        for (int p = 0; p < sys->natoms; ++p) {
                EXPECT_DOUBLE_EQ(sys->r[p].x, 10 * p);
                EXPECT_DOUBLE_EQ(sys->r[p].y, 10 * p);
                EXPECT_DOUBLE_EQ(sys->r[p].z, 10 * p);
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
			
				sys->f = new vec3_t[4]();


                sys->f[proc_id].x = proc_id;
                sys->f[proc_id].y = proc_id;
                sys->f[proc_id].z = proc_id;

                sys->epot = 10.0;
        }

        void TearDown() {
                delete[] sys->f;
			
                delete sys;
        }
};

TEST_F(MPI_reduce_forces_u_test, basic) {

        mpi_reduce_forces(sys);

        if (proc_id == 0) {
                for (int p = 0; p < sys->natoms; ++p) {
                        EXPECT_DOUBLE_EQ(sys->f[p].x, p);
                        EXPECT_DOUBLE_EQ(sys->f[p].y, p);
                        EXPECT_DOUBLE_EQ(sys->f[p].z, p);
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