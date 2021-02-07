#ifndef IO_H
#define IO_H

#include <stdio.h>

#include "physics.h"

#ifdef __cplusplus
extern "C" {
#endif

/* generic file- or pathname buffer length */
#define BLEN 200

// struct to hold the names of all the input files read from stdin
struct _file_names_ {
        char trajfile[BLEN];
        char ergfile[BLEN];
};
typedef struct _file_names_ file_names;

extern int initialise(mdsys_t *sys, FILE *infile, file_names *fnames,
                      int *nprint);
extern void output(mdsys_t *sys, FILE *erg, FILE *traj);

#ifdef __cplusplus
}
#endif

#endif
