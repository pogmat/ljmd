

#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "io.h"
#include "mpi_headers/mpi_test_env.h"

class MPI_brodc_input_params : public ::testing::Test {
      protected:
        int nprocs;
        int proc_id;
        mdsys_t *sys;
        file_names *fnames;

        void SetUp() {
                nprocs = MPITestEnv::get_mpi_procs();
                proc_id = MPITestEnv::get_mpi_rank();
                sys = new mdsys_t;
                fnames = new file_names;
        }

        void TearDown() {
                delete sys;
                delete file_names
        }
};

TEST_F(MPI_brodc_input_params, read_brodc_input_params) {

        FILE *param_file = stdin;

        int nprint;
        char restfile[BLEN];

        read_brodc_input_params(proc_id, sys, param_file, restfile, fnames,
                                nprint);

        if (proc_id == 0) {
                printf("%d\n", sys->natoms);
                printf("%f\n", sys->rcut);
                printf("%s\n", fnames->trajfile);
        }
}

char infile[BLEN];

int main(int argc, char **argv) {

        std::string command_line_arg(argc == 2 ? argv[1] : "");
        testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new MPITestEnv(command_line_arg));

        return RUN_ALL_TESTS();
}
