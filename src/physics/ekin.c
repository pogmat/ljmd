#include "physics.h"
#include <stdio.h>

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_utils.h"
#endif

/* compute kinetic energy */
void ekin(mdsys_t *sys) {
        int i;
        sys->ekin = 0.0;

#if defined(MPI_ENABLED)
        for (i = 0; i < sys->proc_seg->size; ++i) {
#else
        for (i = 0; i < sys->natoms; ++i) {
#endif
                sys->ekin +=
                    0.5 * mvsq2e * sys->mass *
                    (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] +
                     sys->vz[i] * sys->vz[i]);
        }
        sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}
