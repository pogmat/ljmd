#include <stdio.h>
#include <mpi.h>


void mpi_hello( int proc_id) {
		printf("hello from %d\n",proc_id);
}

