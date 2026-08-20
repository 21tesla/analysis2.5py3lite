#ifndef _MYMINIMISATION_H_
#define _MYMINIMISATION_H_

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../inc/myPARAMS.h"
#include "../inc/myDAT.h"
#include "../inc/myRAMACHANDRAN.h"
#include "../inc/mySTRUCT.h"
#include "../inc/myRAMACHANDRAN.h"
#include "../inc/myMACRO.h"

#define D_max 1000              // Max number of dimensions of the search space
#define S_max 500               // Max swarm size
#define R_max 200               // Max number of runs

// Structures

typedef struct {
    gsl_vector     *x;
    double          f;
} position_t;

extern double   PSOMinRotFwd(double (*Func) (const gsl_vector *, void *), gsl_vector * x, void *params);
// extern double   PSOMinRotRev(double (*Func) (const gsl_vector *, void *), gsl_vector * x, void *params, ramaDB_t * ramaDB, polyPeptide_t * polyP, int ppi);
extern double   ConjGradMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params, int iterNb);
extern double   SteepDescMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params);
extern double   SimplexMin(double (*Func_f) (const gsl_vector *, void *), gsl_vector * x, void *params);

#endif
