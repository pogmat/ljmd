#include "gtest/gtest.h"

#include <stdlib.h>

/* inclusion of source file is needed
 * since we want to test a static function  */
extern "C" {
	#include "utils.c"
}

TEST(cleanup_vec3_t_test, empty)
{
	vec3_t *m = NULL;
	cleanup_vec3_t(&m);
	ASSERT_EQ(m, nullptr);
}

TEST(cleanup_vec3_t_test, release)
{
	vec3_t *m = (vec3_t*)malloc(sizeof(vec3_t));
	cleanup_vec3_t(&m);
	ASSERT_EQ(m, nullptr);
}

TEST(cleanup_mdsys_test, initialized)
{
	mdsys_t sys;
	sys.r = (vec3_t *)calloc(2, sizeof(vec3_t));
	sys.v = (vec3_t *)calloc(2, sizeof(vec3_t));
	sys.f = (vec3_t *)calloc(2, sizeof(vec3_t));
	cleanup_mdsys(&sys);
	ASSERT_EQ(sys.r, nullptr);
	ASSERT_EQ(sys.v, nullptr);
	ASSERT_EQ(sys.f, nullptr);
}
