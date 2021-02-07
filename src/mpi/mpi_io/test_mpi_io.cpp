

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

                        sys->rx = new double[4]{0.0, 0.1, 0.2, 0.3};
                        sys->ry = new double[4]{0.0, 0.1, 0.2, 0.3};
                        sys->rz = new double[4]{0.0, 0.1, 0.2, 0.3};

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

        send_pos_vel(nprocs, &proc_seg, sys, vxbuf, vybuf, vzbuf);

        EXPECT_DOUBLE_EQ(sys->rx[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->rx[1], 0.1);
        EXPECT_DOUBLE_EQ(sys->rx[2], 0.2);
        EXPECT_DOUBLE_EQ(sys->rx[3], 0.3);

        EXPECT_DOUBLE_EQ(sys->ry[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->ry[1], 0.1);
        EXPECT_DOUBLE_EQ(sys->ry[2], 0.2);
        EXPECT_DOUBLE_EQ(sys->ry[3], 0.3);

        EXPECT_DOUBLE_EQ(sys->rz[0], 0.0);
        EXPECT_DOUBLE_EQ(sys->rz[1], 0.1);
        EXPECT_DOUBLE_EQ(sys->rz[2], 0.2);
        EXPECT_DOUBLE_EQ(sys->rz[3], 0.3);

        EXPECT_DOUBLE_EQ(sys->vz[0], proc_id);
        EXPECT_DOUBLE_EQ(sys->vy[0], proc_id);
        EXPECT_DOUBLE_EQ(sys->vz[0], proc_id);
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
