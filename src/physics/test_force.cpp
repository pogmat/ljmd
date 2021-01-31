#include "gtest/gtest.h"

/* inclusion of source file is needed
 * since we want to test a static function  */
extern "C" {
	#include "force.c"
}

TEST(pbc_test, inside)
{
	ASSERT_DOUBLE_EQ(pbc(3.0, 10.0 * 0.5), 3.0);
	ASSERT_DOUBLE_EQ(pbc(-2.5, 10.0 * 0.5), -2.5);
}

TEST(pbc_test, outside)
{
	ASSERT_DOUBLE_EQ(pbc(7.0, 10.0 * 0.5), -3.0);
	ASSERT_DOUBLE_EQ(pbc(-8.0, 10.0 * 0.5), 2.0);
}
