#ifndef MPI_UTIL_H
#define MPI_UTIL_H


#include "physics.h"


#ifdef __cplusplus
extern "C"
{
#endif

struct arr_seg {
	int id;
	int idx;
	int size;
	int splitting[3];
};
typedef struct arr_seg arr_seg_t;



void mpi_hello( int proc_id);
	
	
	
extern void init_segments(int nprocs, int proc_id,  arr_seg_t *proc_seg, int size );
extern void split_dimension( const int nprocs,  const  int size , int splitting[3] );

#ifdef __cplusplus	
}
#endif

#endif