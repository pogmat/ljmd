#include <math.h>

#include "physics.h"

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
        double epot = 0.0;

        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);

       #ifdef _OMP_NAIVE
       #pragma omp parallel for default(shared) private(i, j, rx, ry, rz, r, ffac,) reduction(+:epot)
       #endif
        for (i = 0; i < (sys->natoms); ++i) {
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

                                epot += 0.5 * 4.0 * sys->epsilon *
                                             (pow(sys->sigma / r, 12.0) -
                                              pow(sys->sigma / r, 6.0));

                                sys->fx[i] += rx / r * ffac;
                                sys->fy[i] += ry / r * ffac;
                                sys->fz[i] += rz / r * ffac;
                        }
                        sys->epot = epot;
                }
        }
}
