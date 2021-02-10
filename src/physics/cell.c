#include <stdlib.h>

#include "physics.h"

/* (3^3 - 1) / 2 */
#define HALFNEIGH 14

static int pair_number(int ncs) {
	switch (ncs) {
	case 1:
		return 1;
	case 2:
		return 36;
	default:
		return HALFNEIGH * ncs * ncs * ncs;
	}
}

static void build_pairs(mdsys_t *sys) {
	/* Toric fan generating the positive cone */
	int fan[HALFNEIGH][3] = {{ 0,  0,  0},
				 { 1,  0,  0},
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
				 { 1,  1, -1}};

	

	int ncellside = sys->ncellside;
	if (ncellside == 1) {
		sys->cellpairs[0] = 0;
		sys->cellpairs[1] = 0;
		return;
	} 
	if (ncellside == 2) {
		for (int l = 0, p = 0; l < 8; ++l)
			for (int ll = l; ll < 8; ++ll) {
				sys->cellpairs[p++] = l;
				sys->cellpairs[p++] = ll;
			}
		return;
	}
	
	int i, j, k, ii, jj, kk, ll;	
	int *pairptr = sys->cellpairs;
	for (int l = 0; l < ncellside * ncellside * ncellside; ++l) {
		k = l / (ncellside * ncellside);
		j = (l % (ncellside * ncellside)) / ncellside;
		i = (l % (ncellside * ncellside)) % ncellside;
		for (int neig = 0; neig < HALFNEIGH; ++neig) {
			ii = (i + fan[neig][0] + ncellside) % ncellside;
			jj = (j + fan[neig][1] + ncellside) % ncellside;
			kk = (k + fan[neig][2] + ncellside) % ncellside;
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
	int npairs = pair_number(sys->ncellside);
	int ncells = sys->ncellside * sys->ncellside * sys->ncellside;
	sys->cells = (cell_t*)malloc(ncells * sizeof(cell_t));
	sys->cellpairs = (int*)malloc(2 * npairs * sizeof(int));
	build_pairs(sys);
}
