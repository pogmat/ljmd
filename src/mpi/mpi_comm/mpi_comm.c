

#include "io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>

void mpi_send_pos_vel(const int nprocs, const arr_seg_t *proc_seg, mdsys_t *sys,
                      const double *restrict vxbuf,
                      const double *restrict vybuf,
                      const double *restrict vzbuf) {

        /* 	initialise arrays for MPI Scatterv
                count contains the number of elements per segment scattered
                offset indicates the beginning of the segment in the master
           buffer
        */

        int count[nprocs];
        int offsets[nprocs];

        mpi_collective_comm_arrays(nprocs, proc_seg->splitting, count, offsets);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq, vxreq, vyreq, vzreq;

        /*scatter all the positions */

        MPI_Iscatterv(sys->rx, count, offsets, MPI_DOUBLE,
                      &sys->rx[proc_seg->idx], proc_seg->size, MPI_DOUBLE, 0,
                      MPI_COMM_WORLD, &rxreq);
        MPI_Iscatterv(sys->ry, count, offsets, MPI_DOUBLE,
                      &sys->ry[proc_seg->idx], proc_seg->size, MPI_DOUBLE, 0,
                      MPI_COMM_WORLD, &ryreq);
        MPI_Iscatterv(sys->rz, count, offsets, MPI_DOUBLE,
                      &sys->rz[proc_seg->idx], proc_seg->size, MPI_DOUBLE, 0,
                      MPI_COMM_WORLD, &rzreq);

        /* scatter all the velocites*/
        MPI_Iscatterv(vxbuf, count, offsets, MPI_DOUBLE, sys->vx,
                      proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD, &vxreq);
        MPI_Iscatterv(vybuf, count, offsets, MPI_DOUBLE, sys->vy,
                      proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD, &vyreq);
        MPI_Iscatterv(vzbuf, count, offsets, MPI_DOUBLE, sys->vz,
                      proc_seg->size, MPI_DOUBLE, 0, MPI_COMM_WORLD, &vzreq);

        MPI_Wait(&rxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&ryreq, MPI_STATUS_IGNORE);
        MPI_Wait(&rzreq, MPI_STATUS_IGNORE);
        MPI_Wait(&vxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&vyreq, MPI_STATUS_IGNORE);
        MPI_Wait(&vzreq, MPI_STATUS_IGNORE);
}

void mpi_exchange_positions(mdsys_t *sys, const int *count,
                            const int *offsets) {

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq;

        /*	allgather all the positions
                        each process communicates to everyone its own segment
                        and all segments are ressembled into a full positions
           array
                */

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->rx, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &rxreq);

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->ry, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &ryreq);

        MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys->rz, count,
                        offsets, MPI_DOUBLE, MPI_COMM_WORLD, &rzreq);

        MPI_Wait(&rxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&ryreq, MPI_STATUS_IGNORE);
        MPI_Wait(&rzreq, MPI_STATUS_IGNORE);
}

void mpi_reduce_UKT(mdsys_t *sys) {

        MPI_Allreduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &sys->ekin, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &sys->temp, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
}
