

#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

#include "io.h"
#include "mpi_headers/mpi_io.h"
#include "mpi_headers/mpi_test_env.h"

/*
class MPI_brodc_input_params_test : public ::testing::Test {
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
                delete fnames;
        }
};


TEST_F(MPI_brodc_input_params_test, rbip) {

                if (proc_id==0) {
                        freopen("test_params.inp","r",stdin);
                }
        FILE *param_file = stdin;

        int nprint;
        char restfile[BLEN];

                printf("enter param read\n");

        read_brodc_input_params(proc_id, sys, param_file, restfile, fnames,
                                &nprint);

                EXPECT_EQ( sys->natoms , 108 );
                EXPECT_DOUBLE_EQ( sys->mass , 39.948 );
                EXPECT_DOUBLE_EQ( sys->epsilon , 0.2379 );
                EXPECT_DOUBLE_EQ( sys->sigma , 3.405 );
                EXPECT_DOUBLE_EQ( sys->rcut , 8.5 );
                EXPECT_DOUBLE_EQ( sys->box , 17.1580 );
                EXPECT_DOUBLE_EQ( sys->dt , 5.0 );
                EXPECT_EQ( sys->nsteps , 10000 );
                EXPECT_EQ( nprint , 100 );

                EXPECT_TRUE( (strcmp(restfile, "argon_108.rest") == 0) );
                EXPECT_TRUE( (strcmp(fnames->trajfile, "argon_108.xyz") == 0) );
                EXPECT_TRUE( (strcmp(fnames->ergfile, "argon_108.dat") == 0) );



                        printf("%d\n", sys->natoms);
                printf("%f\n", sys->rcut);
                printf("%s\n", fnames->trajfile);

}
*/

int main(int argc, char **argv) {

        for (int i = 0; i < argc; i++) {
                std::cout << i << ":" << argv[i] << std::endl;
        }

        std::string command_line_arg(argc == 2 ? argv[1] : "");
        testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new MPITestEnv(command_line_arg));

        return RUN_ALL_TESTS();
}
