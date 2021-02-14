/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// custom header files
#include "common.h"
#include "io.h"
#include "physics.h"

#if defined(MPI_ENABLED)
#include "mpi_headers/mpi_comm.h"
#include "mpi_headers/mpi_utils.h"
#include <mpi.h>
#endif

/* main */
int main(int argc, char **argv) {
        int nprint;
        file_names fnames;
        FILE *traj, *erg;
        mdsys_t sys;
        double t_start;

#if defined(TIMING)
        double io_t = 0, force_t = 0, e_kin_t = 0, verlet_t = 0, tmp_t;
#endif

#if defined(MPI_ENABLED)
        int nprocs, proc_id, mpi_f;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_f);

        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

        arr_seg_t proc_seg;

        sys.nprocs = nprocs;
        sys.proc_id = proc_id;
        sys.proc_seg = &proc_seg;
        sys.proc_seg->splitting = (int *)malloc(nprocs * sizeof(int));

#endif

#if defined(MPI_ENABLED)
        if (proc_id == 0) {
#endif
                printf("LJMD version %3.1f\n", LJMD_VERSION);
#if defined(MPI_ENABLED)
        }
#endif

        t_start = wallclock();
#if defined(TIMING)
        tmp_t = t_start;
#endif
        /* serial IO */
        char init_err = initialise(&sys, stdin, &fnames, &nprint);
        if (init_err) {
                return init_err;
        }

#if defined(TIMING)
        io_t += wallclock() - tmp_t;
#endif

        /* initialize forces and energies.*/
        sys.nfi = 0;

#if defined(TIMING)
        tmp_t = wallclock();
#endif

#if defined(MPI_ENABLED)
        mpi_broadcast_pos(&sys);
#endif
        force(&sys);
#if defined(MPI_ENABLED)
        mpi_reduce_forces(&sys);

#endif
#if defined(TIMING)
        force_t += wallclock() - tmp_t;
#endif
#if defined(MPI_ENABLED)
        if (proc_id == 0) {
#endif
#if defined(TIMING)
                tmp_t = wallclock();
#endif
                ekin(&sys);
#if defined(TIMING)
                e_kin_t += wallclock() - tmp_t;
#endif

                erg = fopen(fnames.ergfile, "w");
                traj = fopen(fnames.trajfile, "w");

                printf("Startup time: %10.3fs\n", wallclock() - t_start);
                printf("Starting simulation with %d atoms for %d steps.\n",
                       sys.natoms, sys.nsteps);
                printf("     NFI            TEMP            EKIN               "
                       "  EPOT  "
                       "      "
                       "      ETOT\n");
                output(&sys, erg, traj);

#if defined(MPI_ENABLED)
        }
#endif

        /* reset timer */
        t_start = wallclock();

        /**************************************************/
        /* main MD loop */
        for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {

                /* write output, if requested */
#if defined(MPI_ENABLED)
                if (proc_id == 0) {
#endif
                        if ((sys.nfi % nprint) == 0)
                                output(&sys, erg, traj);

                                /* propagate system and recompute energies */
                                /* use the split versin of Verlet algorithm*/
#if defined(TIMING)
                        tmp_t = wallclock();
#endif
                        verlet_1(&sys);
#if defined(TIMING)
                        verlet_t += wallclock() - tmp_t;
#endif
#if defined(MPI_ENABLED)
                }
#endif
#if defined(TIMING)
                tmp_t = wallclock();
#endif
#if defined(MPI_ENABLED)

                mpi_broadcast_pos(&sys);
#endif
                force(&sys);
#if defined(MPI_ENABLED)
                mpi_reduce_forces(&sys);
#endif
#if defined(TIMING)
                force_t += wallclock() - tmp_t;
#endif
#if defined(MPI_ENABLED)
                if (proc_id == 0) {
#endif
#if defined(TIMING)
                        tmp_t = wallclock();
#endif
                        verlet_2(&sys);
#if defined(TIMING)
                        verlet_t += wallclock() - tmp_t;
                        tmp_t = wallclock();
#endif
                        ekin(&sys);
#if defined(TIMING)
                        e_kin_t += wallclock() - tmp_t;
#endif
#if defined(MPI_ENABLED)
                }
#endif
                sleep(0.5);
        }
        /**************************************************/

        /* clean up: close files, free memory */
        cleanup_mdsys(&sys);
#if defined(MPI_ENABLED)
        if (proc_id == 0) {
#endif
                fclose(erg);
                fclose(traj);

#if defined(MPI_ENABLED)
        }
#endif
        /*	print detailed timng information collected
                if MPI we need to reduce the timing values
        */

#if defined(TIMING)
        int t;
        double time_arr[] = {wallclock() - t_start, force_t, verlet_t, e_kin_t};

        for (t = 1; t < 4; ++t) {
                time_arr[t] /= sys.nsteps;
        }

#if defined(MPI_ENABLED)
        double avg_time_arr[] = {0, 0};
        MPI_Reduce(&time_arr, &avg_time_arr, 2, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        if (proc_id == 0) {
                for (t = 0; t < 2; ++t) {
                        avg_time_arr[t] /= nprocs;
                }
#endif

                printf("Simulation Done.\n");
#if defined(MPI_ENABLED)

                printf("Data read time 		: %10.8fs.\n", io_t);
                printf("Avg Force calc time 	: %10.8fs.\n", avg_time_arr[1]);
                printf("Avg Verlet calc time 	: %10.8fs.\n", time_arr[2]);
                printf("Avg Ekin calc time 		: %10.8fs.\n",
                       time_arr[3]);
                printf("Avg total time 		: %10.8fs.\n", avg_time_arr[0]);
        }
#else
        printf("Data read time 		: %10.8fs.\n", io_t);
        printf("Force calc time 	: %10.8fs.\n", time_arr[1]);
        printf("Verlet calc time 	: %10.8fs.\n", time_arr[2]);
        printf("Ekin calc time 		: %10.8fs.\n", time_arr[3]);
        printf("total time 		: %10.3fs.\n", time_arr[0]);
#endif

        /* just print the wallclock time*/
#else
        printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);

#endif

#if defined(MPI_ENABLED)
        MPI_Finalize();
#endif

        return 0;
}
