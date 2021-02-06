#include <stdlib.h>

#include "physics.h"


static void cleanup_vec3_t(vec3_t **p) {
	if(*p) {
		free(*p);
		*p = NULL;
	}
}

void cleanup_mdsys(mdsys_t *sys) {
	cleanup_vec3_t(&sys->r);
	cleanup_vec3_t(&sys->v);
	cleanup_vec3_t(&sys->f);
}
