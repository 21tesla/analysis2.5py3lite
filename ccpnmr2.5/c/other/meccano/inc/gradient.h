// $Id: gradient.h,v 1.1 2008/05/07 10:12:33 wb104 Exp $
#ifndef ALGTOOLS_GRADIENT_H
#define ALGTOOLS_GRADIENT_H 1

// Include files
#include <stdlib.h>
#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"

/** 
 *  Found gradient of function
 *
 *  @author Kirill Miklyaev
 *  @date   2002-05-27
 */

int             gradient(double (*f) (const gsl_vector *, void *), const gsl_vector * v, void *params, gsl_vector * g);

#endif // ALGTOOLS_GRADIENT_H
