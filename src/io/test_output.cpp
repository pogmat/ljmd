#include "io.h"
#include "gtest/gtest.h"

#include <stdio.h>

extern "C" {
#include "output.c"
}

TEST(test_output, output) {

        // Initializing the variables needed for output() function
        mdsys_t sys;
        FILE *erg, *traj;

        // sys = (mdsys_t) {.natoms = 20, .nfi = 101, .ekin = 120.021, .epot =
        // 123.321, .temp = 273.16};

        sys.natoms = 20;
        sys.nfi = 101;
        sys.ekin = 120.021;
        sys.epot = 123.321;
        sys.temp = 273.16;

        sys.rx = (double *)malloc(sys.natoms * sizeof(double));
        sys.ry = (double *)malloc(sys.natoms * sizeof(double));
        sys.rz = (double *)malloc(sys.natoms * sizeof(double));

        for (int i = 0; i < sys.natoms; ++i) {
                sys.rx[i] = (double)i / 100. + 1.0;
                sys.ry[i] = (double)i / 100. + 2.0;
                sys.rz[i] = (double)i / 100. + 3.0;
        }

        erg = fopen("erg_sample.txt", "w");
        traj = fopen("traj_sample.txt", "w");

        output(&sys, erg, traj);

        fclose(erg);
        fclose(traj);

        free(sys.rx);
        free(sys.ry);
        free(sys.rz);

        // Testing the contents of file associated with `erg` filestream
        double total_energy;
        mdsys_t test1_sys;
        erg = fopen("erg_sample.txt", "r");
        fscanf(erg, "%d %lf %lf %lf %lf", &(test1_sys.nfi), &(test1_sys.temp),
               &(test1_sys.ekin), &(test1_sys.epot), &total_energy);
        fclose(erg);

        ASSERT_EQ(test1_sys.nfi, 101);
        ASSERT_DOUBLE_EQ(test1_sys.temp, 273.16);
        ASSERT_DOUBLE_EQ(test1_sys.ekin, 120.021);
        ASSERT_DOUBLE_EQ(test1_sys.epot, 123.321);
        ASSERT_DOUBLE_EQ(total_energy, test1_sys.ekin + test1_sys.epot);
        ASSERT_DOUBLE_EQ(total_energy, sys.ekin + sys.epot);

        // Testing the contents of file associated with `traj` filestream
        mdsys_t test2_sys;
        traj = fopen("traj_sample.txt", "r");
        fscanf(traj, "%d\n nfi=%d etot=%lf\n", &(test2_sys.natoms),
               &(test2_sys.nfi), &total_energy);

        ASSERT_EQ(test2_sys.natoms, 20);
        ASSERT_EQ(test2_sys.nfi, 101);
        ASSERT_DOUBLE_EQ(total_energy, sys.ekin + sys.epot);

        test2_sys.rx = (double *)malloc(test2_sys.natoms * sizeof(double));
        test2_sys.ry = (double *)malloc(test2_sys.natoms * sizeof(double));
        test2_sys.rz = (double *)malloc(test2_sys.natoms * sizeof(double));

        for (int i = 0; i < sys.natoms; ++i) {
                fscanf(traj, "Ar %lf %lf %lf\n", &(test2_sys.rx[i]),
                       &(test2_sys.ry[i]), &(test2_sys.rz[i]));

                ASSERT_DOUBLE_EQ(test2_sys.rx[i], (double)i / 100. + 1.0);
                ASSERT_DOUBLE_EQ(test2_sys.ry[i], (double)i / 100. + 2.0);
                ASSERT_DOUBLE_EQ(test2_sys.rz[i], (double)i / 100. + 3.0);
        }

        fclose(traj);

        free(test2_sys.rx);
        free(test2_sys.ry);
        free(test2_sys.rz);
}
