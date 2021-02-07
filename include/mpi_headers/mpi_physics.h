

#include "mpi_headers/mpi_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void mpi_ekin(arr_seg_t *proc_seg, mdsys_t *sys);
extern void mpi_force(arr_seg_t *proc_seg, mdsys_t *sys);

extern void mpi_verlet_1(arr_seg_t *proc_seg, mdsys_t *sys);
extern void mpi_verlet_2(arr_seg_t *proc_seg, mdsys_t *sys);

#ifdef __cplusplus
}
#endif