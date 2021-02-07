
#include "mpi_headers/mpi_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int mpi_initialise(int nprocs, int proc_id, arr_seg_t *proc_seg, mdsys_t *sys,
                   FILE *infile, file_names *fnames, int *nprint);

int read_brodc_input_params(int proc_id, mdsys_t *sys, FILE *infile,
                            char *restfile, file_names *fnames, int *nprint);

int read_restfile(char *restfile, mdsys_t *sys);

void send_pos_vel(int nprocs, arr_seg_t *proc_seg, mdsys_t *sys, double *vxbuf,
                  double *vybuf, double *vzbuf);

#ifdef __cplusplus
}
#endif