#include <stdlib.h>

#include "physics.h"

/* (3^3 - 1) / 2 */
#define HALFNEIGH 13

/* helper function that determines the number of pairs
 * given a the number fo cells per side.
 * If there are 1 or 2 cell per side every pairs of cells
 * must be considered.
 * If there are 3 or more cells it the fan algorithm applies.
 */
static inline int pair_number(int ncs) {
	switch (ncs) {
	case 1:
		return 0;
	case 2:
		return 28;
	default:
		return HALFNEIGH * ncs * ncs * ncs;
	}
}

/* helper function that build the list of interacting cells */
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
				 { 1,  1, -1}};

	
	/* only one cell per side: the pairs list is trivial */
	int ncellside = sys->ncellside;
	if (ncellside == 1)
		return;
	/* two cells per side: every cell is close to another */
	if (ncellside == 2) {
		for (int l = 0, p = 0; l < 8; ++l)
			for (int ll = l + 1; ll < 8; ++ll) {
				sys->cellpairs[p++] = l;
				sys->cellpairs[p++] = ll;
			}
		return;
	}
	/* more than three cells: the fan algorithm applies */
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

/* this function build the cells */
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

/* helper function that determines in which cell an atom belogs to.
 * Notice that the the cells are set in the following way:
 * 0 => [-box/2                 , -box/2 +     cellsize)
 * 1 => [-box/2 +     cellsize  , -box/2 + 2 * cellsize)
 * 2 => [-box/2 + 2 * cellsize  , -box/2 + 3 * cellsize)
 * ...
 */


inline int which_cell(double x, const mdsys_t *sys) {
	double boxby2 = 0.5 * sys->box;
	return (int)((pbc(x, boxby2) + boxby2) / sys->cellsize);
}
