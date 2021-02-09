#include <stdlib.h>

#include "physics.h"


static void cleanup_double(double **p) {
	if(*p) {
		free(*p);
		*p = NULL;
	}
}

static void cleanup_int(int **p) {
	if (*p) {
		free(*p);
		*p = NULL;
	}
}

static void cleanup_cell_t(cell_t **p) {
	if (*p) {
		free(*p);
		*p = NULL;
	}
}	

void cleanup_mdsys(mdsys_t *sys) {
	cleanup_double(&sys->rx);
	cleanup_double(&sys->ry);
	cleanup_double(&sys->rz);
	cleanup_double(&sys->vx);
	cleanup_double(&sys->vy);
	cleanup_double(&sys->vz);
	cleanup_double(&sys->fx);
	cleanup_double(&sys->fy);
	cleanup_double(&sys->fz);
	cleanup_int(&sys->cellpairs);
	cleanup_cell_t(&sys->cells);
}
