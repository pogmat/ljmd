#include "physics.h"
#include <stdio.h>

/* compute kinetic energy */
void ekin(mdsys_t *sys) {
        int i;
        double ekin = 0.0;

#ifdef _OMP_NAIVE
#pragma omp parallel for default(shared) private(i) reduction(+ : ekin)
#endif
        for (i = 0; i < sys->natoms; ++i) {
                ekin += 0.5 * mvsq2e * sys->mass *
                        (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] +
                         sys->vz[i] * sys->vz[i]);
        }

        sys->ekin = ekin;
        sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}
