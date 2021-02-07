#include <math.h>

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

/* 	force calculations, loop pairwise through the particle position arrays
        ext means that by assumption i belongs to the process
        (i.e. i is in [ proc_seg->idx ; proc_seg->idx + proc_seg->size ] )
        while j is outside this range
*/
static void force_loop_ext(const int i_start, const int i_end,
                           const int j_start, const int j_end, mdsys_t *sys) {
        double r, ffac;
        double rx, ry, rz;
        int i, j;

        for (i = i_start; i < i_end; ++i) {
                for (j = j_start; j < j_end; ++j) {

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

                                sys->fx[i] += rx / r * ffac;
                                sys->fy[i] += ry / r * ffac;
                                sys->fz[i] += rz / r * ffac;
                        }
                }
        }
}

/* 	force calculations, loop pairwise through the particle position arrays
        int means that by assumption both i and j belong to the process
        (i.e. i is in [ proc_seg->idx ; proc_seg->idx + proc_seg->size ] )
*/
static void force_loop_int(const int i_start, const int i_end,
                           const int j_start, const int j_end, mdsys_t *sys) {
        double r, ffac;
        double rx, ry, rz;
        int i, j;

        for (i = i_start; i < i_end; ++i) {
                for (j = i + 1; j < j_end; ++j) {

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

                                sys->epot += 4.0 * sys->epsilon *
                                             (pow(sys->sigma / r, 12.0) -
                                              pow(sys->sigma / r, 6.0));

                                sys->fx[i] += rx / r * ffac;
                                sys->fy[i] += ry / r * ffac;
                                sys->fz[i] += rz / r * ffac;
                                sys->fx[j] += -sys->fx[i];
                                sys->fy[j] += -sys->fy[i];
                                sys->fz[j] += -sys->fz[i];
                        }
                }
        }
}

void mpi_force(arr_seg_t *proc_seg, mdsys_t *sys) {

        /* zero energy and forces */
        sys->epot = 0.0;
        azzero(sys->fx, proc_seg->size);
        azzero(sys->fy, proc_seg->size);
        azzero(sys->fz, proc_seg->size);

        int i_start = proc_seg->idx;
        int i_end = proc_seg->idx + proc_seg->size;

        if (i_start > 0) {
                /*	first run the j loop on the particles
                        with index between 0 and proc_seg->idx
                */
                force_loop_ext(i_start, i_end, 0, proc_seg->idx - 1, sys);
        }

        /*	now run the j loop on the particles
                belonging to the process
        */
        force_loop_int(i_start, i_end, i_start, i_end, sys);

        if (i_end < sys->natoms) {
                /*	finally run the j loop on the particles
                        with index between proc_seg->idx + proc_seg->size and
                   the end
                */
                force_loop_ext(i_start, i_end, proc_seg->idx + 1, sys->natoms,
                               sys);
        }
}
