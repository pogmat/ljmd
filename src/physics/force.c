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
        double epot = 0.0;

        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);
        sys->epot = 0;

        double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
        double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
        double rcsq = sys->rcut * sys->rcut;
        double r6, rinv;
        double r1x, r1y, r1z;
        double f1x, f1y, f1z;

#if defined(MPI_ENABLED)

        for (i = sys->proc_seg->idx;
             i < (sys->proc_seg->idx + sys->proc_seg->size); ++i) {

#else

#ifdef _OMP_NAIVE
#pragma omp parallel for default(shared) private(i, j, rx, ry, rz, ffac, r1x, r1y, r1z, f1x, f1y, f1z, rsq, rinv, r6) reduction(+:epot)
#endif
        for (i = 0; i < (sys->natoms); ++i) {

#endif

                r1x = sys->rx[i];
                r1y = sys->ry[i];
                r1z = sys->rz[i];
                f1x = 0.0;
                f1y = 0.0;
                f1z = 0.0;

#ifdef _OMP_NAIVE
                for (j = 0; j < (sys->natoms); ++j) {
                        if (i == j)
                                continue;
#else
                for (j = i + 1; j < (sys->natoms); ++j) {
#endif

                        /* get distance between particle i and j */
                        rx = pbc(r1x - sys->rx[j], 0.5 * sys->box);
                        ry = pbc(r1y - sys->ry[j], 0.5 * sys->box);
                        rz = pbc(r1z - sys->rz[j], 0.5 * sys->box);
                        rsq = rx * rx + ry * ry + rz * rz;

                        /* compute force and energy if within cutoff */
                        if (rsq < rcsq) {
                                rinv = 1.0 / rsq;
                                r6 = rinv * rinv * rinv;
                                ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
                                epot += r6 * (c12 * r6 - c6);

                                f1x += rx * ffac;
                                f1y += ry * ffac;
                                f1z += rz * ffac;

#ifndef _OMP_NAIVE

#ifdef _OMP_3RD_LAW
#pragma omp critical
                                {
#endif
                                        sys->fx[j] -= rx * ffac;
                                        sys->fy[j] -= ry * ffac;
                                        sys->fz[j] -= rz * ffac;
#ifdef _OMP_3RD_LAW
                                }
#endif

#endif
                        }
                }
                sys->fx[i] += f1x;
                sys->fy[i] += f1y;
                sys->fz[i] += f1z;
        }

        sys->epot = epot;

#ifdef _OMP_NAIVE
        sys->epot *= 0.5;
#endif
}
