#include <stdlib.h>

#include "io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>


void mpi_broadcast_pos(mdsys_t *sys) {

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rreq;

        /*broadcast all the positions */
        MPI_Ibcast(sys->r, 3*sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rreq);

        MPI_Wait(&rreq, MPI_STATUS_IGNORE);
}

void mpi_reduce_forces(mdsys_t *sys) {

        /* reduction in place, two different signatures for master and children
         */
        if (sys->proc_id == 0) {
                MPI_Reduce(MPI_IN_PLACE, sys->f, 3*sys->natoms, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
        } else {
                MPI_Reduce(sys->f, sys->f, 3*sys->natoms, MPI_DOUBLE, MPI_SUM,
                           0, MPI_COMM_WORLD);
        }
        MPI_Allreduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
}
