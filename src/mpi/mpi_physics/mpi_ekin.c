#include "physics.h"
#include <stdio.h>

#include "mpi_headers/mpi_utils.h"

/* compute kinetic energy */
void mpi_ekin(arr_seg_t *proc_seg, mdsys_t *sys) {
        int i;
        sys->ekin = 0.0;

        for (i = 0; i < proc_seg->size; ++i) {
                sys->ekin +=
                    0.5 * mvsq2e * sys->mass *
                    (sys->vx[i] * sys->vx[i] + sys->vy[i] * sys->vy[i] +
                     sys->vz[i] * sys->vz[i]);
        }
        sys->temp = 2.0 * sys->ekin / (3.0 * sys->natoms - 3.0) / kboltz;
}
