#ifndef MPI_UTIL_H
#define MPI_UTIL_H

#include "physics.h"

#ifdef __cplusplus
extern "C" {
#endif

struct arr_seg {
        int id;
        int idx;
        int size;
        int splitting[3];
};
typedef struct arr_seg arr_seg_t;

void mpi_hello(int proc_id);

extern void init_segments(const int nprocs, const int proc_id,
                          arr_seg_t *proc_seg, const int size);
extern void split_dimension(const int nprocs, const int size, int splitting[3]);

extern void mpi_collective_comm_arrays(const int nprocs,
                                       const int *const splitting, int *count,
                                       int *offset);

#ifdef __cplusplus
}
#endif

#endif