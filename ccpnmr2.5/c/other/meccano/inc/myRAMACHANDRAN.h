#ifndef _MYRAMACANDRAN_h_
#define _MYRAMACANDRAN_h_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include "../inc/mySTRUCT.h"
#include "../inc/myGEOMETRY.h"
#include "../inc/mySTRUCT.h"

// typedef struct {
//     int             resolution;
//     gsl_matrix     *gen, *prp, *pro, *gly;
// } ramachandran_t;

typedef struct {
    int             pointNb[NB_AC_AM_TYPE];
    gsl_vector     *phi[NB_AC_AM_TYPE];
    gsl_vector     *psi[NB_AC_AM_TYPE];
} ramaDB_t;

// extern ramachandran_t *InitRama(int resolution);

// extern void     ReadRamaFile(char *fileName, gsl_matrix * ramaMat, int resolution);

// extern void     ReadRamaFiles(ramachandran_t * rama);

// extern double   ERamachandran(ramachandran_t * rama, polyPeptide_t * polyP, int ii);

extern ramaDB_t *ReadRamaDBFile(char *ramaDBFileName);

extern void     AleaRamaDB(ramaDB_t * RamaDB, polyPeptide_t * polyP, int ii, double *phi, double *psi);

extern void     DelRamaDB(ramaDB_t * RamaDB);

#endif
