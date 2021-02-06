#ifndef COMMON_H
#define COMMON_H

#include <sys/time.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* structure to hold 3-vectors */
struct _vec3 {
	double x;
	double y;
	double z;
};
typedef struct _vec3 vec3_t;
	

double wallclock();

void azzero(vec3_t *v, const int n);


#ifdef __cplusplus
}
#endif
	
#endif
