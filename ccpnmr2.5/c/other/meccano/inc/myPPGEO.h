#ifndef _MYPPGEO_H_
#define _MYPPGEO_H_

#include <gsl/gsl_blas.h>
#include "../inc/myGEOMETRY.h"

#define R_CB        1.545
#define R_HA        1.109
#define R_CN        1.334
#define R_CC        1.525
#define R_CH        2.038
#define R_NH        1.020
#define R_CA1H      2.560
#define R_CA2H      2.147

typedef struct {
    gsl_vector     *ca1_ref;
    gsl_vector     *c1__ref;
    gsl_vector     *n2__ref;
    gsl_vector     *hn2_ref;
    gsl_vector     *ca2_ref;
    gsl_vector     *o1__ref;
    gsl_vector     *c_n__ref;
    gsl_vector     *n_ca_ref;
    gsl_vector     *ca_c_ref;
    gsl_vector     *x_pp_ref;
    gsl_vector     *y_pp_ref;
    gsl_vector     *z_pp_ref;
} pPlane_t;

extern pPlane_t *InitPepPla(void);

extern pPlane_t *InitPepPlaFwd(void);

extern pPlane_t *InitPepPlaRev(void);

extern void      DelPepPla(pPlane_t *pPlane);

#endif
