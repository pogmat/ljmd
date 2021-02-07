#include <math.h>
#include <stdio.h>

#include "physics.h"

#include "mpi_headers/mpi_utils.h"

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2) {
        while (x > boxby2)
                x -= 2.0 * boxby2;
        while (x < -boxby2)
                x += 2.0 * boxby2;
        return x;
}

/* compute forces */
void mpi_force(arr_seg_t *proc_seg, mdsys_t *sys) {
        double r, ffac;
        double rx, ry, rz;
        int i, j;

        /* zero energy and forces */
        sys->epot = 0.0;
        azzero(sys->fx, proc_seg->size);
        azzero(sys->fy, proc_seg->size);
        azzero(sys->fz, proc_seg->size);

        int force_idx = 0;

        for (i = proc_seg->idx; i < (proc_seg->idx + proc_seg->size); ++i) {
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
