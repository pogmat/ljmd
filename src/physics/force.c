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
        vec3_t r;
        int i, j;
	
        /* zero energy and forces */
	double epot = 0.0;

        sys->epot = 0.0;
        azzero(sys->f, sys->natoms);


	double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
	double c6  = 4.0 * sys->epsilon * pow(sys->sigma,  6.0);
	double rcsq = sys->rcut * sys->rcut;
	double r6, rinv;

	vec3_t r1;
	vec3_t f1;

	#if defined(MPI_ENABLED)

	for (i = sys->proc_seg->idx;
             i < (sys->proc_seg->idx + sys->proc_seg->size); ++i) {


	#else

        #ifdef _OMP_NAIVE
        #pragma omp parallel for default(shared) private(i, j, r, ffac, r1, f1, rsq, rinv, r6) reduction(+:epot)
        #endif
	for (i = 0; i < (sys->natoms); ++i) {

	#endif

		r1.x = sys->r[i].x;
		r1.y = sys->r[i].y;
		r1.z = sys->r[i].z;
		f1.x = 0.0;
		f1.y = 0.0;
		f1.z = 0.0;


                #ifdef _OMP_NAIVE
                for (j = 0; j < (sys->natoms); ++j) {
                        if(i == j) continue;
                #else

                for (j = i + 1; j < (sys->natoms); ++j) {
                #endif

                        /* get distance between particle i and j */
                        r.x = pbc(sys->r[i].x - sys->r[j].x, 0.5 * sys->box);
                        r.y = pbc(sys->r[i].y - sys->r[j].y, 0.5 * sys->box);
                        r.z = pbc(sys->r[i].z - sys->r[j].z, 0.5 * sys->box);

                        rsq = r.x * r.x + r.y * r.y + r.z * r.z;

			 /* compute force and energy if within cutoff */
                        if (rsq < rcsq) {
				rinv = 1.0 / rsq;
				r6 = rinv * rinv * rinv;
				ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
				epot += r6 * (c12 * r6 - c6);

          
                                f1.x += r.x * ffac;
                                f1.y += r.y * ffac;
                                f1.z += r.z * ffac;

                                #ifndef _OMP_NAIVE

                                #ifdef _OMP_3RD_LAW
                                #pragma omp critical
                                {
                                #endif
        			sys->f[j].x -= r.x * ffac;
				sys->f[j].y -= r.y * ffac;
				sys->f[j].z -= r.z * ffac;
                                #ifdef _OMP_3RD_LAW
                                }
                                #endif

                                #endif

                        }
                        
                }
		sys->f[i].x += f1.x;
		sys->f[i].y += f1.y;
		sys->f[i].z += f1.z;
        }

        sys->epot = epot;

        #ifdef _OMP_NAIVE
        sys->epot *= 0.5;
        #endif
}
