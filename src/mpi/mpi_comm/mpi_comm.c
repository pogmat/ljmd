#include <stdlib.h>

#include "io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>

void mpi_broadcast_pos(mdsys_t *sys) {

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq;

        /*broadcast all the positions */
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rxreq);
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &ryreq);
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rzreq);

        MPI_Wait(&rxreq, MPI_STATUS_IGNORE);
        MPI_Wait(&ryreq, MPI_STATUS_IGNORE);
        MPI_Wait(&rzreq, MPI_STATUS_IGNORE);
}

void mpi_reduce_forces(mdsys_t *sys) {

        /* reduction in place, two different signatures for master and children
         */
        if (sys->proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, sys->fx, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, sys->fy, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, sys->fz, sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
                MPI_Reduce(sys->fx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
                MPI_Reduce(sys->fy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
                MPI_Reduce(sys->fz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
        }
        MPI_Allreduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
}
