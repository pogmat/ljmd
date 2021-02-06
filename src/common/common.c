#include <sys/time.h>
#include "common.h"

/* helper function: get current time in seconds since epoch */

double wallclock() {
        struct timeval t;
        gettimeofday(&t, 0);
        return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

/* helper function: zero out an array */
void azzero(vec3_t *v, const int n) {
        int i;
        for (i = 0; i < n; ++i) {
                (v + i)->x = 0.0;
		(v + i)->y = 0.0;
		(v + i)->z = 0.0;
        }
}
