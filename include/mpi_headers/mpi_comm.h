
#include "physics.h"

#include "mpi_headers/mpi_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void mpi_send_pos_vel(const int nprocs, const arr_seg_t *proc_seg,
                             mdsys_t *sys, const double *vxbuf,
                             const double *vybuf, const double *vzbuf);

extern void mpi_exchange_positions(mdsys_t *sys, const int *count,
                                   const int *offsets);

extern void mpi_reduce_UKT(mdsys_t *sys);

#ifdef __cplusplus
}
#endif