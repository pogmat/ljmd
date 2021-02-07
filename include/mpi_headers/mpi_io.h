
#include "mpi_headers/mpi_utils.h"


int mpi_initialise(int nprocs, int proc_id, arr_seg_t *proc_seg, mdsys_t *sys, FILE *infile, file_names *fnames, int *nprint) ;
int read_brodc_input_params(int nprocs, int proc_id, arr_seg_t *proc_seg, mdsys_t *sys, 
							FILE *infile, char* restfile, file_names *fnames, int *nprint) ;

int read_restfile ( char* restfile, mdsys_t *sys) ;