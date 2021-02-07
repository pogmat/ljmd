#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>

int read_brodc_input_params(int proc_id, mdsys_t *sys, FILE *infile,
                            char *restfile, file_names *fnames, int *nprint) {

        /* prepare two arrays to send the input parameters */

        double send_double_arr[9];
        char send_char_array[2 * BLEN];

        /* Proc 0 reads the input file */
        if (proc_id == 0) {
                char line[BLEN];

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[0] = atoi(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[1] = atof(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[2] = atof(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[3] = atof(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[4] = atof(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[5] = atof(line);

                if (get_a_line(infile, restfile))
                        return 1;
                if (get_a_line(infile, &send_char_array[0]))
                        return 1;
                if (get_a_line(infile, &send_char_array[BLEN]))
                        return 1;
                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[6] = atoi(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[7] = atof(line);

                if (get_a_line(infile, line))
                        return 1;
                send_double_arr[8] = atoi(line);
        }

        /* now broadcast the struct to all processes */
        MPI_Bcast(send_double_arr, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(send_char_array, 2 * BLEN, MPI_CHAR, 0, MPI_COMM_WORLD);

        /* copy all the data into the structs */

        sys->natoms = (int)send_double_arr[0];
        sys->mass = send_double_arr[1];
        sys->epsilon = send_double_arr[2];
        sys->sigma = send_double_arr[3];
        sys->rcut = send_double_arr[4];
        sys->box = send_double_arr[5];
        sys->nsteps = (int)send_double_arr[6];
        sys->dt = send_double_arr[7];
        *nprint = (int)send_double_arr[8];

        memcpy(fnames->trajfile, &send_char_array[0], BLEN * sizeof(char));
        memcpy(fnames->ergfile, &send_char_array[BLEN], BLEN * sizeof(char));
        return 0;
}

int read_restfile(char *restfile, mdsys_t *sys) {

        /* read restart */
        FILE *fp = fopen(restfile, "r");
        int i;

        if (fp) {
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", sys->rx + i, sys->ry + i,
                               sys->rz + i);
                }
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", sys->vx + i, sys->vy + i,
                               sys->vz + i);
                }
                fclose(fp);
                azzero(sys->fx, sys->natoms);
                azzero(sys->fy, sys->natoms);
                azzero(sys->fz, sys->natoms);
        } else {
                return 3;
        }

        return 0;
}

int mpi_initialise(int nprocs, int proc_id, arr_seg_t *proc_seg, mdsys_t *sys,
                   FILE *infile, file_names *fnames, int *nprint) {

        char restfile[BLEN];

        int init_out = read_brodc_input_params(proc_id, sys, infile, restfile,
                                               fnames, nprint);
        if (init_out) {
                perror("Error reading the input file.\n");
                return 1;
        }

        if (proc_id == 1) {
                printf("%d\n", sys->natoms);
                printf("%f\n", sys->rcut);
                printf("%s\n", fnames->trajfile);
        }

        /* initialise the segment structs containing info on this proc's segment
         * of data */
        init_segments(nprocs, proc_id, proc_seg, sys->natoms);

        printf("seg %d idx %d size %d\n", proc_seg->id, proc_seg->idx,
               proc_seg->size);

        /*	allocate the memory
                we doo replicated data and not diustributed, so
                every process requires the full sizes
        */

        sys->rx = (double *)malloc(sys->natoms * sizeof(double));
        sys->ry = (double *)malloc(sys->natoms * sizeof(double));
        sys->rz = (double *)malloc(sys->natoms * sizeof(double));

        sys->fx = (double *)malloc(sys->natoms * sizeof(double));
        sys->fy = (double *)malloc(sys->natoms * sizeof(double));
        sys->fz = (double *)malloc(sys->natoms * sizeof(double));

        /* 	velocity memory is  handled differently than positions or forces
           memory master needs a full-size buffer for input but for simulation
                every process will just need a buffer sized according to number
           of particles
        */

        double *vxbuf = NULL;
        double *vybuf = NULL;
        double *vzbuf = NULL;

        /* first master reads the data
                then we broadcast the positions and scatterv the velocities */
        if (proc_id == 0) {
                sys->vx = (double *)malloc(sys->natoms * sizeof(double));
                sys->vy = (double *)malloc(sys->natoms * sizeof(double));
                sys->vz = (double *)malloc(sys->natoms * sizeof(double));

                int restart_out = read_restfile(restfile, sys);
                if (restart_out) {
                        perror("Cannot read restart file");
                        return 1;
                }

                vxbuf = sys->vx;
                vybuf = sys->vy;
                vzbuf = sys->vz;

        } else {
                sys->vx = NULL;
                sys->vy = NULL;
                sys->vz = NULL;
        }

        sys->vx = (double *)malloc(proc_seg->size * sizeof(double));
        sys->vy = (double *)malloc(proc_seg->size * sizeof(double));
        sys->vz = (double *)malloc(proc_seg->size * sizeof(double));

        /* 	initialise arrays for MPI Scatterv
                count contains the number of elements per segment scattered
                offset indicates the beginning of the segment in the master
           buffer
        */

        int count[nprocs];
        int offsets[nprocs];
        for (int p = 0; p < nprocs; ++p) {
                if (p < proc_seg->splitting[0]) {
                        count[p] = proc_seg->splitting[1];
                        offsets[p] = p * proc_seg->splitting[1];
                } else {
                        count[p] = proc_seg->splitting[2];
                        offsets[p] =
                            proc_seg->splitting[0] * proc_seg->splitting[1] +
                            (p - proc_seg->splitting[0]) *
                                proc_seg->splitting[2];
                }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request rxreq, ryreq, rzreq, vxreq, vyreq, vzreq;

        /*broadcast all the positions */
        MPI_Ibcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rxreq);
        MPI_Ibcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &ryreq);
        MPI_Ibcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD, &rzreq);

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

        if (proc_id == 0) {
                free(vxbuf);
                free(vybuf);
                free(vzbuf);
        }

        return 0;
}