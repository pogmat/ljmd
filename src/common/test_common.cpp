#include "gtest/gtest.h"
#include "common.h"

#include <time.h>


constexpr int SIZE = 20;

class azzero_test: public ::testing::Test {

protected:

	const int size = SIZE;
	double* buf;

	void SetUp()
		{
			buf = new double[size];
			for (int i = 0; i < size; ++i)
				buf[i] = static_cast<double>(i + 1);
		}

	void TearDown()
		{
			delete[] buf;
		}
};


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


TEST_F(azzero_test, doubles)
{
	ASSERT_EQ(size, SIZE);
        for (int i = 0; i < size; ++i)
		ASSERT_DOUBLE_EQ(buf[i], static_cast<double>(i + 1));
	azzero(buf, size);
	for (int i = 0; i < size; ++i)
		ASSERT_DOUBLE_EQ(buf[i], 0.0);
}
