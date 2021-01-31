#include "common.h"
#include "io.h"
#include "gtest/gtest.h"

#include <stdio.h>

extern "C" {
#include "../common/common.c"
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

// TEST(test_input, initialise) {

//   char buffer[] =
//       "108               # natoms\n108               # natoms\n39.948          "
//       "  # mass in AMU\n0.2379            # epsilon in kcal/mol\n3.405         "
//       "    # sigma in angstrom\n8.5               # rcut in angstrom\n17.1580  "
//       "         # box length (in angstrom)\n../src/io/argon_108.rest    # "
//       "restart\nargon_108.xyz     # trajectory\nargon_108.dat     # "
//       "energies\n10000             # nr MD steps\n5.0               # MD time "
//       "step (in fs)\n100               # output print frequency\n";

//   mdsys_t *sys;
//   file_names *fnames;
//   int *nprint;

//   FILE *fstream;

//   fstream = fmemopen(buffer, strlen(buffer), "r");
//   ASSERT_NE(fstream, nullptr);

//   initialise(sys, fstream, fnames, nprint);

//   fclose(fstream);

//   ASSERT_EQ(sys->natoms, 108);
//   ASSERT_DOUBLE_EQ(sys->mass, 39.948);
//   ASSERT_DOUBLE_EQ(sys->epsilon, 0.2379);
//   ASSERT_DOUBLE_EQ(sys->sigma, 3.405);
//   ASSERT_DOUBLE_EQ(sys->rcut, 8.5);
//   ASSERT_DOUBLE_EQ(sys->box, 17.1580);
//   ASSERT_STREQ(fnames->trajfile, "argon_108.xyz");
//   ASSERT_STREQ(fnames->ergfile, "argon_108.dat");
//   ASSERT_EQ(sys->nsteps, 10000);
//   ASSERT_DOUBLE_EQ(sys->dt, 5.0);
//   ASSERT_EQ(*nprint, 100);

//   ASSERT_NE(sys->rx, nullptr);
//   ASSERT_NE(sys->ry, nullptr);
//   ASSERT_NE(sys->rz, nullptr);
//   ASSERT_NE(sys->vx, nullptr);
//   ASSERT_NE(sys->vy, nullptr);
//   ASSERT_NE(sys->vz, nullptr);
//   ASSERT_NE(sys->fx, nullptr);
//   ASSERT_NE(sys->fy, nullptr);
//   ASSERT_NE(sys->fz, nullptr);

//   free(sys->rx);
//   free(sys->ry);
//   free(sys->rz);
//   free(sys->vx);
//   free(sys->vy);
//   free(sys->vz);
//   free(sys->fx);
//   free(sys->fy);
//   free(sys->fz);

// }