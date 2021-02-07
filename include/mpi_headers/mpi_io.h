
#include "mpi_headers/mpi_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int mpi_initialise(const int nprocs, const int proc_id, arr_seg_t *proc_seg,
                   mdsys_t *sys, FILE *infile, file_names *fnames, int *nprint);

int read_brodc_input_params(const int proc_id, mdsys_t *sys, FILE *infile,
                            char *restfile, file_names *fnames, int *nprint);

int read_restfile(const char *restfile, mdsys_t *sys);

#ifdef __cplusplus
}
#endif