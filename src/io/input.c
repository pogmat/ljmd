#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"


/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf) {
        char tmp[BLEN], *ptr;

        /* read a line and cut of comments and blanks */
        if (fgets(tmp, BLEN, fp)) {
                int i;

                ptr = strchr(tmp, '#');
                if (ptr)
                        *ptr = '\0';
                i = strlen(tmp);
                --i;
                while (isspace(tmp[i])) {
                        tmp[i] = '\0';
                        --i;
                }
                ptr = tmp;
                while (isspace(*ptr)) {
                        ++ptr;
                }
                i = strlen(ptr);
                strcpy(buf, tmp);
                return 0;
        } else {
                perror("problem reading input");
                return -1;
        }
        return 0;
}


/* 	handles parsing of initialisation file read through stdin
        allocates the position and velocity arrays and reads
        the restart file with the ICs
*/
int initialise(mdsys_t *sys, FILE *infile ,file_names *fnames, int *nprint) {

        char line[BLEN], restfile[BLEN];
        FILE *fp;
        int i;

        /* read input file */
        if (get_a_line(infile, line))
                return 1;
        sys->natoms = atoi(line);
        if (get_a_line(infile, line))
                return 1;
        sys->mass = atof(line);
        if (get_a_line(infile, line))
                return 1;
        sys->epsilon = atof(line);
        if (get_a_line(infile, line))
                return 1;
        sys->sigma = atof(line);
        if (get_a_line(infile, line))
                return 1;
        sys->rcut = atof(line);
        if (get_a_line(infile, line))
                return 1;
        sys->box = atof(line);
        if (get_a_line(infile, restfile))
                return 1;
        if (get_a_line(infile, fnames->trajfile))
                return 1;
        if (get_a_line(infile, fnames->ergfile))
                return 1;
        if (get_a_line(infile, line))
                return 1;
        sys->nsteps = atoi(line);
        if (get_a_line(infile, line))
                return 1;
        sys->dt = atof(line);
        if (get_a_line(infile, line))
                return 1;
        *nprint = atoi(line);

        /* allocate memory */
	sys->r =  (vec3_t *)malloc(sys->natoms * sizeof(vec3_t));
	sys->v =  (vec3_t *)malloc(sys->natoms * sizeof(vec3_t));
	sys->f =  (vec3_t *)malloc(sys->natoms * sizeof(vec3_t));	

        /* read restart */
        fp = fopen(restfile, "r");
        if (fp) {
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", &(sys->r[i].x), &(sys->r[i].y),
                               &(sys->r[i].z));
                }
                for (i = 0; i < sys->natoms; ++i) {
                        fscanf(fp, "%lf%lf%lf", &(sys->v[i].x), &(sys->v[i].y),
                               &(sys->v[i].z));
                }
                fclose(fp);
                azzero(sys->f, sys->natoms);
        } else {
                perror("cannot read restart file");
                return 3;
        }
        return 0;
}
