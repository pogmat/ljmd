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
        double rsq, ffac;
        double rx, ry, rz;
        int i, j;

        /* zero energy and forces */
        sys->epot = 0.0;
        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);

	double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
	double c6  = 4.0 * sys->epsilon * pow(sys->sigma,  6.0);
	double rcsq = sys->rcut * sys->rcut;
	double r6, rinv;
        for (i = 0; i < (sys->natoms); ++i) {
                for (j = i + 1; j < (sys->natoms); ++j) {

                        /* get distance between particle i and j */
                        rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
                        ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
                        rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
                        rsq = rx * rx + ry * ry + rz * rz;

                        /* compute force and energy if within cutoff */
                        if (rsq < rcsq) {
				rinv = 1.0 / rsq;
				r6 = rinv * rinv * rinv;
				ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
				sys->epot += r6 * (c12 * r6 - c6);

                                sys->fx[i] += rx * ffac;
                                sys->fy[i] += ry * ffac;
                                sys->fz[i] += rz * ffac;
				sys->fx[j] -= rx * ffac;
				sys->fy[j] -= ry * ffac;
				sys->fz[j] -= rz * ffac;
                        }
                }
        }
}
