#include <vector>
#include <utility>
#include <algorithm>

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

class which_cell_test: public ::testing::TestWithParam<std::pair<double, double>> {

protected:
	
	double box = GetParam().first;
	double cellsize = GetParam().second;
	mdsys_t* sys;

	virtual void SetUp() override
		{
		        sys = new mdsys_t;
			sys->box = box;
			sys->cellsize = cellsize;
		}

	virtual void TearDown() override
		{
			delete sys;
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

INSTANTIATE_TEST_SUITE_P(build_pairs_test_parametric,
			 build_pairs_test,
			 ::testing::Values(1, 2, 3, 4, 5));
