#include "gtest/gtest.h"
#include "common.h"

#include <time.h>

TEST(test_wallclock, one_second)
{
	struct timespec t;
	t.tv_sec = 0;
	t.tv_nsec = 100000000;

	double t_test =- wallclock();
	nanosleep(&t, NULL);
	t_test += wallclock();
	ASSERT_TRUE((t_test > 0.09));
	ASSERT_TRUE((t_test < 0.11));
}
