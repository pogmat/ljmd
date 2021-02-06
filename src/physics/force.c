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
        azzero(sys->f, sys->natoms);

	double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
	double c6  = 4.0 * sys->epsilon * pow(sys->sigma,  6.0);
	double rcsq = sys->rcut * sys->rcut;
	double r6, rinv;
        for (i = 0; i < (sys->natoms); ++i) {
                for (j = i + 1; j < (sys->natoms); ++j) {

                        /* get distance between particle i and j */
                        rx = pbc(sys->r[i].x - sys->r[j].x, 0.5 * sys->box);
                        ry = pbc(sys->r[i].y - sys->r[j].y, 0.5 * sys->box);
                        rz = pbc(sys->r[i].z - sys->r[j].z, 0.5 * sys->box);
                        rsq = rx * rx + ry * ry + rz * rz;

                        /* compute force and energy if within cutoff */
                        if (rsq < rcsq) {
				rinv = 1.0 / rsq;
				r6 = rinv * rinv * rinv;
				ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
				sys->epot += r6 * (c12 * r6 - c6);

                                sys->f[i].x += rx * ffac;
                                sys->f[i].y += ry * ffac;
                                sys->f[i].z += rz * ffac;
				sys->f[j].x -= rx * ffac;
				sys->f[j].y -= ry * ffac;
				sys->f[j].z -= rz * ffac;
                        }
                }
        }
}
