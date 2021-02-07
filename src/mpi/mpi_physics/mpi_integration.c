#include "physics.h"

#include "mpi_headers/mpi_utils.h"

// split the verlet algorithm into two functions

/* first part: propagate velocities by half and positions by full step  */
void mpi_verlet_1(arr_seg_t *proc_seg, mdsys_t *sys) {

        int r_idx = proc_seg->idx;

        for (int i = 0; i < proc_seg->size; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
                sys->rx[r_idx] += sys->dt * sys->vx[i];
                sys->ry[r_idx] += sys->dt * sys->vy[i];
                sys->rz[r_idx] += sys->dt * sys->vz[i];
                ++r_idx;
        }
}

/* second part: propagate velocities by another half step */
void mpi_verlet_2(arr_seg_t *proc_seg, mdsys_t *sys) {

        for (int i = 0; i < proc_seg->size; ++i) {
                sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
                sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
                sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        }
}
