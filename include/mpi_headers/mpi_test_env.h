#include "gtest/gtest.h"
#include <mpi.h>

class MPITestEnv : public ::testing::Environment {
      public:
        int arg_nprocs;
        int nprocs;

        explicit MPITestEnv(const std::string &command_line_arg) {
                arg_nprocs = stoi(command_line_arg);
        }

        virtual ~MPITestEnv(){};

      protected:
        virtual void SetUp() {
                char **argv;
                int argc = 0;
		int mpi_f;
                int mpistatus = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_f);
		
                ASSERT_FALSE(mpistatus);

                nprocs = get_mpi_procs();
                ASSERT_EQ(arg_nprocs, nprocs);

                int proc_id;
                MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
        }

        virtual void TearDown() {
                int mpistatus = MPI_Finalize();
                ASSERT_FALSE(mpistatus);
        }

      public:
        static int get_mpi_procs() {
                int procs;
                MPI_Comm_size(MPI_COMM_WORLD, &procs);
                return procs;
        }

        static int get_mpi_rank() {
                int id;
                MPI_Comm_rank(MPI_COMM_WORLD, &id);
                return id;
        }
};
