#include <math.h>

#include "physics.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_utils.h"
#endif

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2) {
        while (x > boxby2)
                x -= 2.0 * boxby2;
        while (x < -boxby2)
                x += 2.0 * boxby2;
        return x;
}

/* compute forces */
void force(mdsys_t *sys) {
        double r, ffac;
        double rx, ry, rz;
        int i, j;

        /* zero energy and forces */
        sys->epot = 0.0;

        int force_idx = 0;
#if defined(MPI_ENABLED)
        azzero(sys->fx, sys->proc_seg->size);
        azzero(sys->fy, sys->proc_seg->size);
        azzero(sys->fz, sys->proc_seg->size);

        for (i = sys->proc_seg->idx;
             i < (sys->proc_seg->idx + sys->proc_seg->size); ++i) {
#else
        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);
        for (i = 0; i < (sys->natoms); ++i) {
#endif
                for (j = 0; j < (sys->natoms); ++j) {

                        /* particles have no interactions with themselves */
                        if (i == j)
                                continue;

                        /* get distance between particle i and j */
                        rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
                        ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
                        rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
                        r = sqrt(rx * rx + ry * ry + rz * rz);

                        /* compute force and energy if within cutoff */
                        if (r < sys->rcut) {
                                ffac = -4.0 * sys->epsilon *
                                       (-12.0 * pow(sys->sigma / r, 12.0) / r +
                                        6 * pow(sys->sigma / r, 6.0) / r);

                                sys->epot += 0.5 * 4.0 * sys->epsilon *
                                             (pow(sys->sigma / r, 12.0) -
                                              pow(sys->sigma / r, 6.0));

                                sys->fx[force_idx] += rx / r * ffac;
                                sys->fy[force_idx] += ry / r * ffac;
                                sys->fz[force_idx] += rz / r * ffac;
                        }
                }
                ++force_idx;
        }
}
