#include "common.h"
#include "io.h"
#include "gtest/gtest.h"

#include <stdio.h>

extern "C" {
#include "common.c"
#include "input.c"
}

TEST(test_input, get_a_line) {
  char line[BLEN], buffer[BLEN];
  FILE *fstream;

  // Testing integer parsing
  sprintf(buffer, "42          # This is an integer input\n");
  fstream = fmemopen(buffer, BLEN, "r");
  ASSERT_NE(fstream, nullptr);

  get_a_line(fstream, line);
  ASSERT_EQ(42, atoi(line));
  fclose(fstream);

  // Testing double precision number parsing
  sprintf(buffer, "3.14159265359# This is an integer input\n");
  fstream = fmemopen(buffer, BLEN, "r");
  ASSERT_NE(fstream, nullptr);

  get_a_line(fstream, line);
  ASSERT_DOUBLE_EQ(3.14159265359, atof(line));
  fclose(fstream);

  // Testing string parsing
  sprintf(buffer, "this_is_a_file_name.txt     # This is an integer input\n");
  fstream = fmemopen(buffer, BLEN, "r");
  ASSERT_NE(fstream, nullptr);

  get_a_line(fstream, line);
  ASSERT_STREQ("this_is_a_file_name.txt", line);
  fclose(fstream);
}

TEST(test_input, initialise) {

  char buffer[] ="108               # natoms\n"
  "39.948            # mass in AMU\n"
  "0.2379            # epsilon in kcal/mol\n"
  "3.405             # sigma in angstrom\n"
  "8.5               # rcut in angstrom\n"
  "17.1580           # box length (in angstrom)\n"
  "test_template.txt    # restart\n"
  "argon_108.xyz     # trajectory\n"
  "argon_108.dat     # energies\n"
  "10000             # nr MD steps\n"
  "5.0               # MD time step (in fs)\n"
  "100               # output print frequency\n";

  int nprint;
  file_names fnames;
  mdsys_t sys;

  FILE *fstream;

  fstream = fmemopen(buffer, strlen(buffer), "r");
  ASSERT_NE(fstream, nullptr);

  FILE* template_file;
  template_file = fopen("test_template.txt", "w");
  for(int i=0; i<2*108; ++i){
    fprintf(template_file, "%f\t%f\t%f\n", (double)i/100., (double)i/100., (double)i/100.);
  }
  fclose(template_file);

  int return_val = initialise(&sys, fstream, &fnames, &nprint);

  fclose(fstream);

  ASSERT_EQ(sys.natoms, 108);
  ASSERT_DOUBLE_EQ(sys.mass, 39.948);
  ASSERT_DOUBLE_EQ(sys.epsilon, 0.2379);
  ASSERT_DOUBLE_EQ(sys.sigma, 3.405);
  ASSERT_DOUBLE_EQ(sys.rcut, 8.5);
  ASSERT_DOUBLE_EQ(sys.box, 17.1580);
  ASSERT_STREQ(fnames.trajfile, "argon_108.xyz");
  ASSERT_STREQ(fnames.ergfile, "argon_108.dat");
  ASSERT_EQ(sys.nsteps, 10000);
  ASSERT_DOUBLE_EQ(sys.dt, 5.0);
  ASSERT_EQ(nprint, 100);

  for(int i=0; i<108; ++i){
    ASSERT_DOUBLE_EQ(sys.rx[i], (double)i/100.);
    ASSERT_DOUBLE_EQ(sys.ry[i], (double)i/100.);
    ASSERT_DOUBLE_EQ(sys.rz[i], (double)i/100.);

    ASSERT_DOUBLE_EQ(sys.vx[i], (double)(i+108)/100.);
    ASSERT_DOUBLE_EQ(sys.vy[i], (double)(i+108)/100.);
    ASSERT_DOUBLE_EQ(sys.vz[i], (double)(i+108)/100.);
  }
  
  ASSERT_NE(sys.rx, nullptr);
  ASSERT_NE(sys.ry, nullptr);
  ASSERT_NE(sys.rz, nullptr);
  ASSERT_NE(sys.vx, nullptr);
  ASSERT_NE(sys.vy, nullptr);
  ASSERT_NE(sys.vz, nullptr);
  ASSERT_NE(sys.fx, nullptr);
  ASSERT_NE(sys.fy, nullptr);
  ASSERT_NE(sys.fz, nullptr);

  free(sys.rx);
  free(sys.ry);
  free(sys.rz);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);

}
