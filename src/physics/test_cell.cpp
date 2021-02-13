#include <vector>
#include <utility>
#include <algorithm>
#include <random>
#include <map>

#include <iostream>

#include "gtest/gtest.h"

extern "C" {
	#include "cell.c"
}

inline int min(int x, int y) {return x < y ? x : y;}
inline int max(int x, int y) {return x < y ? y : x;}

static int pbc_int(int x, int size) {
	while (x > size / 2)
		x -= size;
	while (x < - size / 2)
		x += size;
	return x;
}

class build_pairs_test: public ::testing::TestWithParam<int> {

protected:
	
	mdsys_t sys;
	int size = GetParam();
	int size3 = size * size * size;
	
	virtual void SetUp() override
		{
			sys.ncellside = size;
			sys.cellpairs = new int[2 * pair_number(size)];
		}

	virtual void TearDown() override
		{
			delete[] sys.cellpairs;
		}
};


class which_cell_test: public ::testing::TestWithParam<std::pair<int, double>> {

protected:

	mdsys_t sys;
	int ncell = GetParam().first;
	double boxsize = GetParam().second;
	double cellsize = boxsize / ncell;
	using real_dist = std::uniform_real_distribution<double>;
	std::map<int, real_dist> dists;
	std::default_random_engine re;

	virtual void SetUp() override
		{
			sys.box = boxsize;
			sys.cellsize = cellsize;

			for (int i = -2 * ncell; i < 2 * ncell; ++i)
				dists.emplace(i, real_dist(-boxsize / 2 + i * cellsize,
						           -boxsize / 2 + (i + 1) * cellsize));   
			
		}
	
};

TEST_P(build_pairs_test, fill)
{
	using vec_t = std::vector<std::pair<int, int>>;
	vec_t there, here;
	build_pairs(&sys);

	int x, y, M, m, npairs;

	npairs = pair_number(size);
	
	for (int l = 0; l < 2 * npairs; l += 2) {
		x = sys.cellpairs[l];
		y = sys.cellpairs[l + 1];
		M = max(x, y);
		m = min(x, y);
		there.emplace_back(m, M);
	}
	
	std::sort(there.begin(), there.end());

	int i, j, k, ii, jj, kk, di, dj, dk;
	for (int l = 0; l < size3; ++l) {
		k = l / (size * size);
		j = (l % (size * size)) / size;
		i = (l % (size * size)) % size;
		for (int ll = l + 1; ll < size3; ++ll) {
			kk = ll / (size * size);
			jj = (ll % (size * size)) / size;
			ii = (ll % (size * size)) % size;
			di = pbc_int(i - ii, size);
			dj = pbc_int(j - jj, size);
			dk = pbc_int(k - kk, size);
			if (di * di + dj * dj + dk * dk < 4) {
				here.emplace_back(l, ll);
			}
		}
	}
}

TEST_P(which_cell_test, normal) {
	int idx;
	real_dist d;
	//std::cout << boxsize << std::endl;
	//std::cout << ncell << std::endl;
	//std::cout << cellsize << std::endl;
	for (const auto& el: dists) {
	        idx = el.first;
		d = el.second;
		for (int l = 0; l < 5; ++l)
			ASSERT_EQ(which_cell(d(re), &sys), (idx + 4 * ncell) % ncell);
	}
}


INSTANTIATE_TEST_SUITE_P(build_pairs_test_parametric,
			 build_pairs_test,
			 ::testing::Values(1, 2, 3, 4, 5));


INSTANTIATE_TEST_SUITE_P(which_cell_test_parametric,
			 which_cell_test,
			 ::testing::Values(std::make_pair<int, double>(1, 10.0),
					   std::make_pair<int, double>(2, 16.0),
					   std::make_pair<int, double>(3, 5.64),
					   std::make_pair<int, double>(4, 23.5),
					   std::make_pair<int, double>(5, 4.33)));
