#include "gtest/gtest.h"

#include <stdlib.h>

/* inclusion of source file is needed
 * since we want to test a static function  */
extern "C" {
	#include "utils.c"
}

TEST(cleanup_double_test, empty)
{
	double *m = NULL;
	cleanup_double(&m);
	ASSERT_EQ(m, nullptr);
}

TEST(cleanup_double_test, release)
{
	double *m = (double*)malloc(sizeof(double));
	cleanup_double(&m);
	ASSERT_EQ(m, nullptr);
}

TEST(cleanup_mdsys_test, initialized)
{
	mdsys_t sys;
	sys.rx = (double*)calloc(2, sizeof(double));
	sys.ry = (double*)calloc(2, sizeof(double));
	sys.rz = (double*)calloc(2, sizeof(double));
	sys.vx = (double*)calloc(2, sizeof(double));
	sys.vy = (double*)calloc(2, sizeof(double));
	sys.vz = (double*)calloc(2, sizeof(double));
	sys.fx = (double*)calloc(2, sizeof(double));
	sys.fy = (double*)calloc(2, sizeof(double));
	sys.fz = (double*)calloc(2, sizeof(double));
	cleanup_mdsys(&sys);
	ASSERT_EQ(sys.rx, nullptr);
	ASSERT_EQ(sys.ry, nullptr);
	ASSERT_EQ(sys.rz, nullptr);
	ASSERT_EQ(sys.vx, nullptr);
	ASSERT_EQ(sys.vy, nullptr);
	ASSERT_EQ(sys.vz, nullptr);
	ASSERT_EQ(sys.fx, nullptr);
	ASSERT_EQ(sys.fy, nullptr);
	ASSERT_EQ(sys.fz, nullptr);
}
