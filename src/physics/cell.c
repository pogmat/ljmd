#include <stdlib.h>

#include "physics.h"

/* (3^3 - 1) / 2 */
#define HALFNEIGH 13


static void build_pairs(mdsys_t *sys) {
	/* Toric fan generating the positive cone */
	int fan[HALFNEIGH][3] = {{ 1,  0,  0},
				{ 0,  1,  0},
				{ 0,  0,  1},
				{ 1,  1,  0},
				{ 1,  0,  1},
				{ 0,  1,  1},
				{ 1,  1,  1},
				{ 1, -1,  0},
				{ 1,  0, -1},
				{ 0,  1, -1},
				{-1,  1,  1},
				{ 1, -1,  1},
				{ 1,  1,  1}};

	int i, j, k, ii, jj, kk, ll;
	int ncellside = sys->ncellside;
	int *pairptr = sys->cellpairs;
	for (int l = 0; l < ncellside * ncellside * ncellside; ++l) {
		k = l / (ncellside * ncellside);
		j = (l % (ncellside * ncellside)) / ncellside;
		i = (l % (ncellside * ncellside)) % ncellside;
		for (int neig = 0; neig < HALFNEIGH; ++neig) {
			ii = (i + fan[neig][0]) % ncellside;
			jj = (j + fan[neig][1]) % ncellside;
			kk = (k + fan[neig][2]) % ncellside;
			ll = ii + jj * ncellside + kk * ncellside * ncellside; 
			*(pairptr++) = l;
			*(pairptr++) = ll;
		}
	}
	
}


void build_cells(mdsys_t *sys) {
	const double tolerance = 1.01;
	
	sys->ncellside = (int)(sys->box / (sys->rcut * tolerance));
	sys->cellsize = sys->box / sys->ncellside;
	int ncell = sys->ncellside * sys->ncellside * sys->ncellside;
	sys->cells = (cell_t*)malloc(ncell * sizeof(cell_t));
	sys->cellpairs = (int*)malloc(2 * HALFNEIGH * ncell * sizeof(int));
	build_pairs(sys);
}
