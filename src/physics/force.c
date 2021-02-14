#include <math.h>
#include <stdio.h>

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
        double rsq, ffac;
        double rx, ry, rz;
        int i, j;

        /* zero energy and forces */
        sys->epot = 0.0;
        double pot_energy;
        azzero(sys->f, sys->natoms);

        double force_x, force_y, force_z;

#if defined(MPI_ENABLED)
        for (i = sys->proc_seg->idx;
             i < (sys->proc_seg->idx + sys->proc_seg->size); ++i) {
                for (j = i + 1; j < (sys->natoms); ++j) {
#else
        for (i = 0; i < (sys->natoms); ++i) {
                for (j = 0; j < (sys->natoms); ++j) {

                        /* particles have no interactions with themselves */
                        if (i == j)
                                continue;

#endif

                        /* get distance between particle i and j */
                        rx = pbc(sys->r[i].x - sys->r[j].x, 0.5 * sys->box);
                        ry = pbc(sys->r[i].y - sys->r[j].y, 0.5 * sys->box);
                        rz = pbc(sys->r[i].z - sys->r[j].z, 0.5 * sys->box);
                        rsq = sqrt(rx * rx + ry * ry + rz * rz);

                        /* compute force and energy if within cutoff */
                        if (rsq < sys->rcut) {
                                ffac =
                                    -4.0 * sys->epsilon *
                                    (-12.0 * pow(sys->sigma / rsq, 12.0) / rsq +
                                     6 * pow(sys->sigma / rsq, 6.0) / rsq);

                                pot_energy = 4.0 * sys->epsilon *
                                             (pow(sys->sigma / rsq, 12.0) -
                                              pow(sys->sigma / rsq, 6.0));

#if !defined(MPI_ENABLED)
                                pot_energy /= 2;
#endif
                                sys->epot += pot_energy;

                                force_x = rx / rsq * ffac;
                                force_y = ry / rsq * ffac;
                                force_z = rz / rsq * ffac;

                                sys->f[i].x += force_x;
                                sys->f[i].y += force_y;
                                sys->f[i].z += force_z;
#if defined(MPI_ENABLED)
                                sys->f[j].x -= force_x;
                                sys->f[j].y -= force_y;
                                sys->f[j].z -= force_z;
#endif
                        }
                }
        }
}
