#include "gtest/gtest.h"
#include <math.h>

#include "mpi_headers/mpi_utils.h"


TEST(test_split_dimension, all)
{
	int size = 1000;
	int nprocs = GetParam();
	
	int splitting[3];
	
	split_dimension( nprocs, size, splitting);
	
	float base = ((float)size)/((float)nprocs);
	int h1 = (int)ceil(base);
	int h2 = (int)floor(base);
	
	
	ASSERT( splitting[0] > 0);	
	EXPECT( splitting[0] <= nprocs);
	
	EXPECT( splitting[1] <= h1);
	EXPECT( splitting[1] >= h2);
	EXPECT( splitting[2] <= h1);
	EXPECT( splitting[3] >= h2);
	
	EXPECT_EQ( splitting[0]*splitting[1] + (nprocs - splitting[0])*splitting[2] , size);
					
}



INSTANTIATE_TEST_SUITE_P(test_split_dimension_p,
			 test_split_dimension,
			 ::testing::Values(0, 1, 2, 3));