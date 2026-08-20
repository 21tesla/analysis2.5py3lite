// $Id: gradient.c,v 1.4 2010-09-09 16:04:50 wb104 Exp $
// Include files

// local
#include "../inc/gradient.h"

//-----------------------------------------------------------------------------
// Implementation file
// 
// 2002-05-26 : Kirill Miklyaev
//-----------------------------------------------------------------------------

int gradient(double (*f) (const gsl_vector *, void *), const gsl_vector * v, void *params, gsl_vector * g) {
    unsigned int    i, j, k;
    double          h = GSL_SQRT_DBL_EPSILON;
    double          a[4], d[4], a3;
    double          result, abserr;
    static double   x;

    gsl_vector     *v1 = gsl_vector_alloc(v->size);

    gsl_vector_memcpy(v1, v);

    for (j = 0; j < v->size; ++j) {
        for (i = 0; i < 4; i++) {
            x = gsl_vector_get(v1, j);
            a[i] = x + (i - 2.0) * h;

            gsl_vector_set(v1, j, a[i]);
            d[i] = f(v1, params);
            gsl_vector_set(v1, j, x);

        }

        for (k = 1; k < 5; k++) {
            for (i = 0; i < 4 - k; i++) {
                d[i] = (d[i + 1] - d[i]) / (a[i + k] - a[i]);
            }
        }

        a3 = fabs(d[0] + d[1] + d[2] + d[3]);

        if (a3 < 100.0 * GSL_SQRT_DBL_EPSILON) {
            a3 = 100.0 * GSL_SQRT_DBL_EPSILON;
        }

        h = pow(GSL_SQRT_DBL_EPSILON / (2.0 * a3), 1.0 / 3.0);

        if (h > 100.0 * GSL_SQRT_DBL_EPSILON) {
            h = 100.0 * GSL_SQRT_DBL_EPSILON;
        }

        x = gsl_vector_get(v1, j);

        gsl_vector_set(v1, j, x + h);
        double          f1 = f(v1, params);

        gsl_vector_set(v1, j, x);

        gsl_vector_set(v1, j, x - h);
        double          f2 = f(v1, params);

        gsl_vector_set(v1, j, x);

        result = (f1 - f2) / (2.0 * h);
        abserr = fabs(100.0 * a3 * h * h);

        gsl_vector_set(g, j, result);

    }

    gsl_vector_free(v1);

    return GSL_SUCCESS;

}
