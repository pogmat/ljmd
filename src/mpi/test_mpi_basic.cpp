

#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "mpi_headers/mpi_test_env.h"

class MPITest : public ::testing::Test {
      protected:
        int nprocs;

        void SetUp() { nprocs = MPITestEnv::get_mpi_procs(); }
};

TEST_F(MPITest, hello) {

        int this_rank = MPITestEnv::get_mpi_rank();

        if (this_rank == 0) {
                std::cout << "from " << this_rank << " we have " << nprocs
                          << " procs running" << std::endl;
        }
}

TEST_F(MPITest, reduce) {

        int this_rank = MPITestEnv::get_mpi_rank();
        int sendbuf = 1;
        int recvbuf = 0;

        MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (this_rank == 0) {
                ASSERT_EQ(nprocs, recvbuf);
                std::cout << "result of reduction : " << recvbuf << std::endl;
        } else {
                ASSERT_EQ(0, recvbuf);
        }
}

int main(int argc, char **argv) {

        for (int i = 0; i < argc; i++) {
                std::cout << i << ":" << argv[i] << std::endl;
        }

        std::string command_line_arg(argc == 2 ? argv[1] : "");
        testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new MPITestEnv(command_line_arg));

        return RUN_ALL_TESTS();
}
