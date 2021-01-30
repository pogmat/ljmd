#ifndef COMMON_H
#define COMMON_H

#include <sys/time.h>

#ifdef __cplusplus
extern "C"
{
#endif


double wallclock();

void azzero(double *d, const int n);


#ifdef __cplusplus
}
#endif
	
#endif
