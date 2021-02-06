#include "physics.h"

/* compute kinetic energy */
void ekin(mdsys_t *sys) {
        int i;

        sys->ekin = 0.0;
        for (i = 0; i < sys->natoms; ++i) {
                sys->ekin +=
                    0.5 * mvsq2e * sys->mass *
                    (sys->v[i].x * sys->v[i].x + sys->v[i].y * sys->v[i].y +
                     sys->v[i].z * sys->v[i].z);
        }
        sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}
