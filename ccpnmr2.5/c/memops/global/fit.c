/*
======================COPYRIGHT/LICENSE START==========================

fit.c: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

*/
#include "fit.h"

#include "cpmg.h"
#include "line_fit.h"
#include "nonlinear_model.h"
#include "random.h"

#define  MAX_MODEL_ITER  20
#define  MAX_CONDITION  4
/*
#define  MAX_CONDITION  10
*/
#define  CHISQ_STOP_CRITERION  (1.0e-1)
 
#define  MIN_DELTA_X  (1.0e-2)
#define  SMALL_X  (1.0e-7)

#define  SMALL_ETA  (1.0e-7)

/*
#define  BSCALE  10.0
*/
#define  BSCALE  1.0

#define  CHECK_NONLINEAR_MODEL(stage) \
	 {   if (nonlinear_model(x, y, fit->w, n, params_fit, fit->covar, \
				fit->alpha, fit->beta, fit->da, fit->ap, fit->dy_da, \
				fit->piv, fit->row, fit->col, nparams, &chisq, \
				&lambda, f, stage, NULL, error_msg) == CCPN_ERROR) \
	         RETURN_ERROR_MSG("singular data");   }

typedef struct Log_linear_fit
{
    FIT_STRUCT
    int nalloc;
    int *used_indices;
    float *xx;
    float *yy;
}   *Log_linear_fit;

typedef struct Non_linear_fit
{
    FIT_STRUCT
    int nalloc;
    float *w;
    float **covar;
    float **alpha;
    float *beta;
    float *da;
    float *ap;
    float *dy_da;
    int *piv;
    int *row;
    int *col;
}   *Non_linear_fit;

static void free_log_linear_memory(Log_linear_fit fit)
{
    FREE(fit->used_indices, int);
    FREE(fit->xx, float);
    FREE(fit->yy, float);

    fit->nalloc = 0;
}

static CcpnStatus alloc_log_linear_memory(Log_linear_fit fit, int n, CcpnString error_msg)
{
    sprintf(error_msg, "allocating log linear memory");
 
    if (n < fit->nalloc)
	return  CCPN_OK;
 
    if (fit->nalloc > 0)
    {
	free_log_linear_memory(fit);
    }
    else
    {
	fit->used_indices = NULL;
	fit->xx = NULL;
	fit->yy = NULL;
    }

    MALLOC(fit->used_indices, int, n);
    MALLOC(fit->xx, float, n);
    MALLOC(fit->yy, float, n);
 
    fit->nalloc = n;

    return  CCPN_OK;
}

static void free_nonlinear_memory(Non_linear_fit fit)
{
    static int nparams = 3;

    FREE(fit->w, float);

    FREE(fit->beta, float);
    FREE(fit->da, float);
    FREE(fit->ap, float);
    FREE(fit->dy_da, float);
    FREE2(fit->covar, float, nparams);
    FREE2(fit->alpha, float, nparams);

    FREE(fit->piv, int);
    FREE(fit->row, int);
    FREE(fit->col, int);
}

static CcpnStatus alloc_nonlinear_memory(Non_linear_fit fit, int n, CcpnString error_msg)
{
    static int nparams = 4;
    int i;

    sprintf(error_msg, "allocating nonlinear memory");
 
    if (n < fit->nalloc)
	return  CCPN_OK;
 
    if (fit->nalloc > 0)
	free_nonlinear_memory(fit);

    MALLOC(fit->w, float, n);
    for (i = 0; i < n; i++)
	fit->w[i] = 1; /* weight */

    MALLOC(fit->beta, float, nparams);
    MALLOC(fit->da, float, nparams);
    MALLOC(fit->ap, float, nparams);
    MALLOC(fit->dy_da, float, nparams);
    MALLOC2(fit->covar, float, nparams, nparams);
    MALLOC2(fit->alpha, float, nparams, nparams);

    MALLOC(fit->piv, int, nparams);
    MALLOC(fit->row, int, nparams);
    MALLOC(fit->col, int, nparams);
 
    fit->nalloc = n;

    return  CCPN_OK;
}

static CcpnStatus fit_linear(int n, float *x, float *y, float noise,
		float *params_fit, float *params_dev,
		float *y_fit, float *chisq, CcpnString error_msg)
{
    float *sigma = NULL, a, b, std_a, std_b, corr_ab, goodness;

/*
int i;
printf("n = %d\n", n);
for (i = 0; i < n; i++) printf("  %d: x = %4.3e y = %4.3e\n", i, x[i], y[i]);
*/
    CHECK_STATUS(line_fit(n, x, y, sigma, y_fit, &a, &b, &std_a, &std_b, &corr_ab, &goodness, error_msg));
/*
for (i = 0; i < n; i++) printf("  %d: yfit = %4.3e\n", i, y_fit[i]);
*/

    params_fit[0] = b; /* NOTE a and b backwards! */
    params_fit[1] = a;

/* since passed sigma=NULL to line_fit, goodness is just chi2 of fit */
    if (n > 2)
      *chisq = goodness / ((n-2) * noise * noise);
    else
      *chisq = 0;

    if (params_dev)
    {
	params_dev[0] = std_b; /* NOTE a and b backwards! */
	params_dev[1] = std_a;
    }

    return  CCPN_OK;
}

static CcpnStatus fit_log_linear(Log_linear_fit fit, int n, float *x, float *y,
			float noise, float *params_fit, float *params_dev,
			float *y_fit, float *chisq, CcpnString error_msg)
{
    int i, j, imax, nused;
    float ymax, s, *sigma = NULL, a, b, std_a, std_b, corr_ab, goodness, d, d2;
    
    CHECK_STATUS(alloc_log_linear_memory(fit, n, error_msg));

    imax = 0;
    ymax = ABS(y[0]);

    for (i = 1; i < n; i++)
    {
	if (ABS(y[i]) > ymax)
	{
	    imax = i;
	    ymax = ABS(y[i]);
	}
    }

    if (ymax == 0)
	RETURN_ERROR_MSG("y's all zero");

    if (y[imax] > 0)
	s = 1;
    else
	s = -1;

    nused = 0;
    for (i = j = 0; i < n; i++)
    {
	fit->yy[j] = s * y[i];

	if (fit->yy[j] > 0)
	{
	    fit->yy[j] = (float) log((double) fit->yy[j]);
	    fit->xx[j] = x[i];
	    fit->used_indices[j] = i;
	    j++;
	    nused++;
	}
    }

    if (nused < 2)
	RETURN_ERROR_MSG("not enough valid data points");

    CHECK_STATUS(line_fit(nused, fit->xx, fit->yy, sigma, y_fit, &a, &b, &std_a, &std_b, &corr_ab, &goodness, error_msg));

    params_fit[0] = s * (float) exp((double) a);
    params_fit[1] = - b;

    for (i = 0; i < nused; i++)
    {
	y_fit[i] = (float) exp((double) y_fit[i]);
	j = fit->used_indices[i];
	d = y[j] - y_fit[i];
	d2 = d * d;
    }

    if (nused > 2)
      *chisq = d2 / ((nused-2) * noise * noise);
    else
      *chisq = 0;

    if (params_dev)
    {
	/* first one is approximate */
	params_dev[0] = std_a * (float) exp((double) 2*a);
	params_dev[1] = std_b;
    }

    if (nused < n) /* put back in missing terms, zeroed (arbitrary) */
    {
	for (i = nused-1; i >= 0; i--)
	{
	    y_fit[fit->used_indices[i]] = y_fit[i];

	    if (i > 0)
	    {
		for (j = fit->used_indices[i-1]+1; j < fit->used_indices[i]; j++)
		    y_fit[j] = 0;
	    }
	}
    }

    return  CCPN_OK;
}

static void nonlinear_func(float x, int nparams, float *params, float *y, float *dy_da)
{
    float a = params[0], b = params[1], c = 0;
    float s = (float) exp((double) (-b*x));
    float t = a * s;

    if (nparams == 3)
    {
	c = params[2];
	dy_da[2] = 1;
    }

    *y = t + c;
    dy_da[0] = s;
    dy_da[1] = - x * t;
}

static void nonlinear2_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    nonlinear_func(x, 2, params, y, dy_da);
}

static void nonlinear3_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    nonlinear_func(x, 3, params, y, dy_da);
}

static void slow_exchange_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1], c = params[2];
    double s, t;

/*
    printf("slow_exchange_func: a = %4.3e, b = %4.3e, c = %4.3e\n", a, b, c);
*/
    dy_da[2] = 1;

    x *= b;
    x = ABS(x);

    if (x < SMALL_X)
    {
	s = x * x / 6;
	*y = (float) (a*s + c);
	dy_da[0] = (float) s;
	dy_da[1] = (float) (2 * a * s / b);
    }
    else
    {
	s = sin((double) x) / x;
	t = cos((double) x);
	*y = (float) (a*(1-s) + c);
	dy_da[0] = (float) (1 - s);
	dy_da[1] = (float) a * (s - t) / b;
    }
}

static void langmuir_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], s;

/*
    printf("langmuir_func: a = %4.3e\n", a);
*/
    x *= a;
    s = 1 + x;

    if (s != 0)
    {
	*y = x / s;
	dy_da[0] = *y / (a*s);
    }
    else
    {
	*y = LARGE_FLOAT;
	dy_da[0] = LARGE_FLOAT;
    }
}

static void kd_shift_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
/* was below until 8 Aug 2007; re-instated as an alternative 2 Jun 2008
    float a = params[0], b = params[1], s, t;

    t = b + x;
    s = t*t - 4*x;
    if (s > 0)
	s = sqrt(s);
    else
	s = 0; // arbitrary
    dy_da[0] = t - s;
    *y = a * dy_da[0];
    if (s != 0)
	dy_da[1] = a * (1 - t/s);
    else
	dy_da[1] = LARGE_FLOAT;
*/
    float a = params[0], b = params[1], c = params[2], s, t, w;

    if (x > 0)
    {
	t = 1 + b/(4*x);
	s = t*t - 1;
	if (s > 0)
	{
	    s = sqrt(s);
	    w =  t - s - c;
	    *y =  a * w;
	    dy_da[0] = w;
	    dy_da[1] = a * (1 - t/s) / (4*x);
	}
	else
	{
	    *y = dy_da[0] = dy_da[1] = 0; /* arbitrary */
	}

	dy_da[2] = -a;
    }
    else
    {
	*y = dy_da[0] = dy_da[1] = 0; /* arbitrary */
	dy_da[2] = -a;
    }
}

static void kd_alt_shift_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1], s, t;

    t = b + x;
    s = t*t - 4*x;
    if (s > 0)
	s = sqrt(s);
    else
	s = 0; // arbitrary
    dy_da[0] = t - s;
    if (s != 0)
    {
        *y = a * dy_da[0];
	dy_da[1] = a * (1 - t/s);
    }
    else
    {
	*y = LARGE_FLOAT;
	dy_da[1] = -LARGE_FLOAT;
    }
}

static void inversion_recovery_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1], s, t;

    s = (float) exp((double) (- b*x/BSCALE));
    t = 0.5 - s;
    *y = a * t;
    dy_da[0] = t;
    dy_da[1] = a * x * s / BSCALE;
}

static void gaussian_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1];
    float x2 = x * x;
    float s = (float) exp((double) (-b*x2));
    float t = a * s;

    *y = t;
    dy_da[0] = s;
    dy_da[1] = - x2 * t;
}

static void cosine_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1];
    float s = (float) cos((double) (b*x));
    float t = a * s;
    float u = (float) sin((double) (b*x));

    *y = t;
    dy_da[0] = s;
    dy_da[1] = - a * x * u;
}

static void cpmg3_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
/*
    see Mulder, Mittermaier, Hon, Dahlquist and Kay
    Nature Structural Biology 8 (2001) 932-935
    http://www.nature.com/nsmb/journal/v8/n11/pdf/nsb1101-932.pdf
    x = nu = 1 / (2 tau)
*/
    float R2max = params[0], kex = params[1], dw = params[2];
    float psi, psis, Dp, Dm, etap, etam, v, vv;
    float cosh_etap, sinh_etap, cos_etam, sin_etam;
    float dD_dkex, dv_dkex, dD_ddw, dv_ddw;
    float detap_dkex, detam_dkex, detap_ddw, detam_ddw;

    psi = kex*kex - dw*dw;
    if (psi > 0)
    {
	Dp = kex*kex / psi;
        Dm = dw*dw / psi;
	psis = sqrt(psi);
        etap = psis / (2.0*x);
        cosh_etap = cosh(etap);
        sinh_etap = sinh(etap);
	v = Dp*cosh_etap - Dm;
        vv = sqrt(v*v - 1);

	dD_dkex = -2*kex*Dm / psi;
	detap_dkex = kex / (2.0*x*psis);
	dv_dkex = dD_dkex * (cosh_etap - 1) + Dp * sinh_etap * detap_dkex;

	dD_ddw = 2*dw*Dm / psi;
	detap_ddw = -dw / (2.0*x*psis);
	dv_ddw = dD_ddw * (cosh_etap - 1) + Dp * sinh_etap * detap_ddw;
    }
    else
    {
	psi = -psi;
        Dp = dw*dw / psi;
	Dm = kex*kex / psi;
	psis = sqrt(psi);
        etam = psis / (2.0*x);
        cos_etam = cos(etam);
        sin_etam = sin(etam);
	v = Dp - Dm*cos_etam;
        vv = sqrt(v*v - 1);

	dD_dkex = 2*kex*Dp / psi;
	detam_dkex = -kex / (2.0*x*psis);
	dv_dkex = dD_dkex * (1 - cos_etam) - Dm * sin_etam * detam_dkex;

	dD_ddw = -2*dw*Dm / psi;
	detam_ddw = dw / (2.0*x*psis);
	dv_ddw = dD_ddw * (1 - cos_etam) - Dm * sin_etam * detam_ddw;
    }

    *y = R2max + 0.5*kex - x*acosh(v);
    dy_da[0] = 1;
    dy_da[1] = 0.5 - x*dv_dkex/vv;
    dy_da[2] = - x*dv_ddw/vv;
}

static void cpmg4_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
/*
    see Mulder, Mittermaier, Hon, Dahlquist and Kay
    Nature Structural Biology 8 (2001) 932-935
    http://www.nature.com/nsmb/journal/v8/n11/pdf/nsb1101-932.pdf
    x = nu = 1 / (2 tau)
*/
    float R2max = params[0], kAB = params[1], kBA = params[2], dw = params[3];
    float psi, zeta, Dp, Dm, etap, etam, v, vv, t, s, kex;
    float cosh_etap, sinh_etap, cos_etam, sin_etam;
    float dD_dkAB, dv_dkAB, dD_dkBA, dv_dkBA, dD_ddw, dv_ddw, fkAB, fkBA, fdw, gp, gm, hh;
    float detap_dkAB, detam_dkAB, detap_dkBA, detam_dkBA, detap_ddw, detam_ddw;
    float params3[3];

    kex = kAB + kBA;
    psi = kex*kex - dw*dw;
    zeta = 2*dw*(kAB - kBA);
    t = sqrt(psi*psi + zeta*zeta);
    s = (psi + 2*dw*dw) / t;

    Dp = 0.5 * (1 + s);
    Dm = 0.5 * (-1 + s);
    etap = sqrt(0.5*(psi+t)) / (2.0*x);
    etam = sqrt(0.5*(-psi+t)) / (2.0*x);

    if ((ABS(etap) < SMALL_ETA) || (ABS(etam) < SMALL_ETA))
    {
	params3[0] = R2max;
	params3[1] = kAB + kBA;
	params3[2] = dw;
        cpmg3_func(x, params3, y, dy_da, user_data);
	dy_da[3] = dy_da[2];
	dy_da[2] = dy_da[1];

	return;
    }

    cosh_etap = cosh(etap);
    sinh_etap = sinh(etap);
    cos_etam = cos(etam);
    sin_etam = sin(etam);
    v = Dp*cosh_etap - Dm*cos_etam;
    vv = sqrt(v*v - 1);
    fkAB = (psi*kex + zeta*dw) / t;
    fkBA = (psi*kex - zeta*dw) / t;
    fdw = (-psi*dw + zeta*(kAB-kBA)) / t;
    gp = 1.0 / (8*x*x*etap);
    gm = 1.0 / (8*x*x*etam);
    hh = cosh_etap - cos_etam;

    dD_dkAB = kex/t - s*fkAB/t;
    dD_dkBA = kex/t - s*fkBA/t;
    dD_ddw = dw/t - s*fdw/t;
    detap_dkAB = gp * (kex + fkAB);
    detap_dkBA = gp * (kex + fkBA);
    detap_ddw = gp * (-dw + fdw);
    detam_dkAB = gm * (-kex + fkAB);
    detam_dkBA = gm * (-kex + fkBA);
    detam_ddw = gm * (dw + fdw);
    dv_dkAB = hh * dD_dkAB + Dp * sinh_etap * detap_dkAB + Dm * sin_etam * detam_dkAB;
    dv_dkBA = hh * dD_dkBA + Dp * sinh_etap * detap_dkBA + Dm * sin_etam * detam_dkBA;
    dv_ddw =  hh * dD_ddw  + Dp * sinh_etap * detap_ddw  + Dm * sin_etam * detam_ddw;

    *y = R2max + 0.5*kex - x*acosh(v);
    dy_da[0] = 1;
    dy_da[1] = 0.5 - x*dv_dkAB/vv;
    dy_da[2] = 0.5 - x*dv_dkBA/vv;
    dy_da[3] = - x*dv_ddw/vv;
}

static void nonlinear_inverse_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1], c = params[2];
    float s = (float) exp((double) (-b*x));
    float t = a * s;

    *y = c - t;
    dy_da[0] = - s;
    dy_da[1] = x * t;
    dy_da[2] = 1;
}

static void nonlinear_inverse2_func(float x, float *params, float *y, float *dy_da, void *user_data)
{
    float a = params[0], b = params[1];
    float s = (float) exp((double) (-b*x));
    float t = 1.0 - s;

    *y = a * t;
    dy_da[0] = t;
    dy_da[1] = x * a * s;
}

/* finds smallest non-negative value of x and y at corresponding x */
static void find_y_xmin(float *x, float *y, int n, float *xmin, float *y_at_xmin)
{
    int i;

    *xmin = x[0];
    *y_at_xmin = y[0];

    for (i = 0; i < n; i++)
    {
	if (x[i] > 0)
	{
	    *xmin = x[i];
	    *y_at_xmin = y[i];
	    return;
	}
    }

    for (i = 0; i < n; i++)
    {
	if ((x[i] > 0) && (x[i] < *xmin))
	{
	    *xmin = x[i];
	    *y_at_xmin = y[i];
	}
    }
}

static void init_params_fit(int method, int n,
			float *x, float *y, float *params_fit)
{
    int i, n1, n2;
    float y1, y2, abs_y;
    float xmin, xmax, y_at_xmin, y_at_xmax, ymax, x_at_ymax, ymin;
    float xmax2, y_at_xmax2;
    Bool done;

    if (method == SLOW_EXCHANGE_FIT)
    {
/*
  let (a, b, c) = params_fit
  we are fitting y = c + a (1 - sin(bx)/(bx))
  limit as x --> 0 of y is c
  limit as x --> infinity of y is c+a
  maximum of y is first bx > 0 where tan(bx) = bx, i.e. bx ~ 4.5
  estimate c = (y where x is minimum)
  estimate a = (y where x is maximum) - c
  estimate b = 4.5 / (x where y is maximum)
*/
	xmin = xmax = x_at_ymax = x[0];
	y_at_xmin = y_at_xmax = ymax = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] < xmin)
	    {
		xmin = x[i];
		y_at_xmin = y[i];
	    }
	    else if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }

	    if (y[i] > ymax)
	    {
		ymax = y[i];
		x_at_ymax = x[i];
	    }
	}

	params_fit[0] = y_at_xmax - y_at_xmin;
	if (x_at_ymax > 0) /* should always be true for slow exchange */
	    params_fit[1] = 4.5 / x_at_ymax;
	else
	    params_fit[1] = 1;
	params_fit[2] = y_at_xmin;
	return;
    }
    else if (method == LANGMUIR_FIT)
    {
/*
  let (a,) = params_fit
  we are fitting y = ax / (1 + ax)
  estimate a = y/x for minimum x
*/
	find_y_xmin(x, y, n, &xmin, &y_at_xmin);

	if (xmin > 0)
	    params_fit[0] = y_at_xmin / xmin;
	else
	    params_fit[0] = 1; /* arbitrary */

	return;
    }
    else if (method == KD_SHIFT_FIT)
    {
/* was this until 8 Aug 2007; re-instated as an alternative 2 Jun 2008
	xmax = x[0];
	for (i = 1; i < n; i++)
	    xmax = MAX(xmax, x[i]);

	params_fit[0] = 0.5 * xmax;
	params_fit[1] = 1.1;
*/
	xmax = x[0];
	y_at_xmax = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }
	}

	params_fit[0] = y_at_xmax;
	params_fit[1] = 0.5 * xmax; /* not sure if this is sensible */
	params_fit[2] = 0;
	return;
    }
    else if (method == KD_ALT_SHIFT_FIT)
    {
	xmax = x[0];
	y_at_xmax = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }
	}

/***
	xmin = xmax;
	y_at_xmin = y_at_xmax;
	for (i = 0; i < n; i++)
	{
	    if (x[i] > 0 && x[i] < xmin)
	    {
		xmin = x[i];
		y_at_xmin = y[i];
	    }
	}
***/

/*
	params_fit[0] = 0.5 * xmax;
*/

/**
	done = CCPN_FALSE;
	if (xmax > 1)
	{
	    xmax2 = -LARGE_FLOAT;
	    for (i = 0; i < n; i++)
	    {
		if ((x[i] > xmax2) && (x[i] != xmax))
		{
		    xmax2 = x[i];
		    y_at_xmax2 = y[i];
		    done = CCPN_TRUE;
		}
	    }
	}
**/

/**
	if (done && (xmax2 > 1))
	{
	    params_fit[1] = (y_at_xmax-y_at_xmax2)/(y_at_xmax2/xmax2 - y_at_xmax/xmax);
	    params_fit[0] = 0.5*y_at_xmax*(1+xmax/params_fit[1]);
	}
	else
**/
	{
	    params_fit[0] = y_at_xmax;
/***
	    params_fit[0] = 0.5 * y_at_xmax;
	    if (y_at_xmin > 0 && xmin < 0.2)
	    {
		params_fit[1] = 2 * params_fit[0] * xmin / y_at_xmin;
	    }
	    else
***/
	    {
	    	params_fit[1] = 1.1;
	    }
	}

	return;
    }
    else if (method == INVERSION_RECOVERY_FIT)
    {
	xmax = x[0];
	y_at_xmax = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }
	}

	params_fit[0] = 2 * y_at_xmax;
	params_fit[1] = 0;
	return;
    }
    else if (method == GAUSSIAN_FIT)
    {
	xmin = xmax = x[0];
	y_at_xmax = y_at_xmin = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }
	    else if (x[i] < xmin)
	    {
		xmin = x[i];
		y_at_xmin = y[i];
	    }
	}

	params_fit[1] = (float) log((double) ABS(y_at_xmax/y_at_xmin)) / (xmax*xmax - xmin*xmin);
	params_fit[0] = y_at_xmin / ((float) exp((double) (-params_fit[1]*xmin)));
	return;
    }
    else if (method == COSINE_FIT)
    {
	xmax = ABS(x[0]);
	ymax = y[0];
	for (i = 1; i < n; i++)
	{
	    xmax = MAX(xmax, ABS(x[i]));
	    ymax = MAX(ymax, y[i]);
	}

	params_fit[0] = ymax;
        if (xmax > 0)
	    params_fit[1] = 1.0 / xmax;
        else
	    params_fit[1] = 1.0;
	return;
    }
    else if (method == CPMG3_FAST_FIT)
    { 
	cpmg3_fast_init_params(n, x, y, params_fit);
/*
	printf("R2max init = %.2f\n", params_fit[0]);
 	printf("kex init = %.2f\n", params_fit[1]);
	printf("dw init = %.2f\n", params_fit[2]);
*/

	return;
    }
    else if (method == CPMG3_SLOW_FIT)
    { 
	cpmg3_slow_init_params(n, x, y, params_fit);
/*
	printf("R2max init = %.2f\n", params_fit[0]);
 	printf("kex init = %.2f\n", params_fit[1]);
	printf("dw init = %.2f\n", params_fit[2]);
*/

	return;
    }
    else if (method == CPMG4_FAST_FIT)
    { 
	cpmg4_fast_init_params(n, x, y, params_fit);
/*
	printf("R2max init = %.2f\n", params_fit[0]);
 	printf("kAB init = %.2f\n", params_fit[1]);
 	printf("kBA init = %.2f\n", params_fit[2]);
	printf("dw init = %.2f\n", params_fit[3]);
*/

	return;
    }
    else if (method == CPMG4_SLOW_FIT)
    { 
	cpmg4_slow_init_params(n, x, y, params_fit);
/*
	printf("R2max init = %.2f\n", params_fit[0]);
 	printf("kAB init = %.2f\n", params_fit[1]);
 	printf("kBA init = %.2f\n", params_fit[2]);
	printf("dw init = %.2f\n", params_fit[3]);
*/

	return;
    }
    else if (method == NONLINEAR_INVERSE_FIT)
    {
        n1 = n2 = -1;
        y1 = y2 = -1;
        ymax = 0;

        for (i = 0; i < n; i++)
        {
	    abs_y = ABS(y[i]);
            if (abs_y > ABS(ymax))
                ymax = y[i];

	    if (n1 >= 0)
	    {
	        if (ABS(x[i] - x[n1]) < MIN_DELTA_X)
		    continue;
	    }

	    if (n2 >= 0)
	    {
	        if (ABS(x[i] - x[n2]) < MIN_DELTA_X)
		    continue;
	    }

	    if (abs_y < y1)
	    {
	        n2 = n1;
	        y2 = y1;
	        n1 = i;
	        y1 = abs_y;
	    }
	    else if (abs_y < y2)
	    {
	        n2 = i;
	        y2 = abs_y;
	    }
        }

	y1 = ymax - y1;
	y2 = ymax - y2;

        if (n2 < 0)
	{
	    params_fit[0] = y1;
	    params_fit[1] = 0;
	    params_fit[2] = ymax;
	}
	else
	{
	    params_fit[1] = (float) log((double) ABS(y1/y2)) / (x[n2] - x[n1]);
	    params_fit[0] = y[n1] / ((float) exp((double) (-params_fit[1]*x[n1])));
	    params_fit[2] = ymax;
	}

	return;
    }
    else if (method == NONLINEAR_INVERSE2_FIT)
    {
	xmin = xmax = x[0];
	y_at_xmin = y_at_xmax = y[0];
	for (i = 1; i < n; i++)
	{
	    if (x[i] < xmin)
	    {
		xmin = x[i];
		y_at_xmin = y[i];
	    }
	    else if (x[i] > xmax)
	    {
		xmax = x[i];
		y_at_xmax = y[i];
	    }
	}

	params_fit[0] = y_at_xmax;
        if ((xmin > 0) && (y_at_xmax != 0) && (ABS(y_at_xmin) < (0.25 * ABS(y_at_xmax))))
	    params_fit[1] = y_at_xmin / (y_at_xmax * xmin);
	else
	    params_fit[1] = 0.0;

	return;
    }

    n1 = n2 = -1;
    y1 = y2 = -1;

    for (i = 0; i < n; i++)
    {
	if (n1 >= 0)
	{
	    if (ABS(x[i] - x[n1]) < MIN_DELTA_X)
		continue;
	}

	if (n2 >= 0)
	{
	    if (ABS(x[i] - x[n2]) < MIN_DELTA_X)
		continue;
	}

	abs_y = ABS(y[i]);

	if (abs_y > y1)
	{
	    n2 = n1;
	    y2 = y1;
	    n1 = i;
	    y1 = abs_y;
	}
	else if (abs_y > y2)
	{
	    n2 = i;
	    y2 = abs_y;
	}
    }

    if (n2 < 0)
    {
	params_fit[0] = y1;
	params_fit[1] = 0;
    }
    else
    {
	/* in theory y[n1] and y[n2] could have opposite signs */
	/* but very unlikely so ignore this */
	params_fit[1] = (float) log((double) ABS(y1/y2)) / (x[n2] - x[n1]);
	params_fit[0] = y[n1] / ((float) exp((double) (-params_fit[1]*x[n1])));
    }

    if (method >= NONLINEAR3_FIT)
	params_fit[2] = 0;

/*
    printf("params_fit = %4.3e, %4.3e, %4.3e\n", params_fit[0], params_fit[1], params_fit[2]);
*/
}

static float calculate_fit(int method, float x, float *params_fit)
{
    float a, b, c, s, t;

    a = params_fit[0];

    if (method == LANGMUIR_FIT)
    {
	x *= a;
	if (x != -1)
	    return x / (1 + x);
	else
	    return LARGE_FLOAT; /* arbitrary */
    }

    b = params_fit[1];

    if (method == KD_SHIFT_FIT)
    {
/* was this until 8 Aug 2007; re-instated as an alternative 2 Jun 2008
	t = b + x;
	s = t*t - 4*x;
	if (s > 0)
	    s = sqrt(s);
	else
	    s = 0; // arbitrary
	return a * (t-s);
*/
	if (x > 0)
	{
	    c = params_fit[2];
	    t = 1 + b/(4*x);
	    s = t*t - 1;
	    if (s > 0)
	        s = sqrt(s);
	    else
	        s = 0; /* arbitrary */

	    return a * ((t-s) - c);
	}
	else
	{
	    return 0; /* arbitrary */
	}
    }

    if (method == KD_ALT_SHIFT_FIT)
    {
	t = b + x;
	s = t*t - 4*x;
	if (s > 0)
	    s = sqrt(s);
	else
	    s = 0; // arbitrary
/*
printf("fit: x=%.3f, a=%.3f, b=%.3f, s=%.3f, y=%.3f\n", x, a, b, s, a*(t-s));
*/
	return a * (t-s);
    }

    if (method == INVERSION_RECOVERY_FIT)
	return ((float) a * (0.5 - exp((double) (-b*x/BSCALE))));

    if (method == GAUSSIAN_FIT)
	return ((float) a * exp((double) (-b*x*x)));

    if (method == COSINE_FIT)
	return ((float) a * cos((double) (b*x)));

    if ((method == CPMG3_FAST_FIT) || (method == CPMG3_SLOW_FIT))
	return cpmg3(x, params_fit);

    if ((method == CPMG4_FAST_FIT) || (method == CPMG4_SLOW_FIT))
	return cpmg4(x, params_fit);

    if (method == NONLINEAR_INVERSE_FIT)
    {
	c =  params_fit[2];
	return (c - (float) a * exp((double) (-b*x)));
    }

    if (method == NONLINEAR_INVERSE2_FIT)
	return a * (1.0 - (float) exp((double) (-b*x)));

    if (method >= NONLINEAR3_FIT)
	c =  params_fit[2];
    else
	c = 0;

    if (method == LINEAR_FIT)
    {
	return a*x + b;
    }
    else if (method <= NONLINEAR3_FIT)
    {
	return ((float) a * exp((double) (-b*x))) + c;
    }
    else if (method == SLOW_EXCHANGE_FIT)
    {
	x *= b;
	x = ABS(x);
	if (x < SMALL_X)
	    return c + a*x*x/6;
	else
	    return (float) a * (1 - sin((double) (x))/x) + c;
    }

    return 0.0; /* should never be here but keep compiler happy */
}

static CcpnStatus fit_nonlinear(Non_linear_fit fit, int n, float *x, float *y,
	 		float noise, float *params_fit, float *params_dev,
			float *y_fit, float *p_chisq, CcpnString error_msg)
{
    int i, cond, iter, nparams = get_method_nparams(fit->method);
    float lambda, chisq, old_chisq, chisq_stop_criterion;/*, ysum;*//*, chisq_scale;*/
    Bool singular;
    Nonlinear_model_func f;

/*
printf("method = %s, n = %d\n", fit_methods[method], n);
*/

    CHECK_STATUS(alloc_nonlinear_memory(fit, n, error_msg));

/*
    params_fit[0] = y[0];
    params_fit[1] = 0;
*/
    init_params_fit(fit->method, n, x, y, params_fit);

    if (fit->method == NONLINEAR2_FIT)
    {
	f = nonlinear2_func;
    }
    else if (fit->method == NONLINEAR3_FIT)
    {
	f = nonlinear3_func;
/*
	params_fit[2] = 0;
*/
    }
    else if (fit->method == SLOW_EXCHANGE_FIT)
    {
	f = slow_exchange_func;
    }
    else if (fit->method == LANGMUIR_FIT)
    {
	f = langmuir_func;
    }
    else if (fit->method == KD_SHIFT_FIT)
    {
	f = kd_shift_func;
    }
    else if (fit->method == INVERSION_RECOVERY_FIT)
    {
	f = inversion_recovery_func;
    }
    else if (fit->method == KD_ALT_SHIFT_FIT)
    {
	f = kd_alt_shift_func;
    }
    else if (fit->method == GAUSSIAN_FIT)
    {
	f = gaussian_func;
    }
    else if (fit->method == COSINE_FIT)
    {
	f = cosine_func;
    }
    else if ((fit->method == CPMG3_FAST_FIT) || (fit->method == CPMG3_SLOW_FIT))
    {
	f = cpmg3_func;
    }
    else if ((fit->method == CPMG4_FAST_FIT) || (fit->method == CPMG4_SLOW_FIT))
    {
	f = cpmg4_func;
    }
    else if (fit->method == NONLINEAR_INVERSE_FIT)
    {
	f = nonlinear_inverse_func;
    }
    else if (fit->method == NONLINEAR_INVERSE2_FIT)
    {
	f = nonlinear_inverse2_func;
    }

/*
    chisq_scale = 100.0 / n;
    chisq_stop_criterion = chisq_scale * CHISQ_STOP_CRITERION;
*/
    chisq_stop_criterion = CHISQ_STOP_CRITERION * noise * noise;
/*
    printf("chisq_stop_criterion = %4.3e\n", chisq_stop_criterion);
*/
 
/*
    ysum = 0;
    for (i = 0; i < n; i++)
	ysum += MAX(ABS(y[i]), noise);
    if (ysum > 0)
    {
        for (i = 0; i < n; i++)
	    fit->w[i] = ABS(y[i]) / ysum;
    }
*/

    CHECK_NONLINEAR_MODEL(INITIAL_STAGE);
 
    for (iter = cond = 0; (iter < MAX_MODEL_ITER) &&
				(cond < MAX_CONDITION); iter++)
    {
/*
	printf("iteration %d\n", iter);
*/
 
	old_chisq = chisq;
 
	CHECK_NONLINEAR_MODEL(GENERAL_STAGE);
 
	if (chisq > old_chisq)
	    cond = 0;
	else if ((old_chisq - chisq) < chisq_stop_criterion)
	    cond++;

/*
printf("iter = %d, cond = %d, lambda = %4.3e, chisq = %4.3e, old_chisq = %4.3e, params0 = %4.3e, params1 = %4.3e",
                iter, cond, lambda, chisq, old_chisq, params_fit[0], params_fit[1]);
if (fit->method >= NONLINEAR3_FIT)
printf(", params2 = %4.3e", params_fit[2]);
printf("\n");
*/
    }

/*
printf("chisq = %4.3e, params0 = %4.3e, params1 = %4.3e\n",
                chisq, params_fit[0], params_fit[1]);
for (i = 0; i < n; i++)
        printf("i = %d, x = %4.3e, y = %4.3e, yfit = %4.3e\n", i, x[i], y[i],
                (float) (params_fit[0]*exp((double) (-params_fit[1]*x[i]))));
int j;
for (i = 0; i < nparams; i++)
for (j = 0; j < nparams; j++)
{
printf("covar[%d][%d] = %e\n", i, j, fit->covar[i][j]);
printf("alpha[%d][%d] = %e\n", i, j, fit->alpha[i][j]);
}
*/

    if (iter == MAX_MODEL_ITER)
	RETURN_ERROR_MSG("fit did not converge");

    if (params_dev)
	CHECK_NONLINEAR_MODEL(FINAL_STAGE);

    for (i = 0; i < n; i++)
	y_fit[i] = calculate_fit(fit->method, x[i], params_fit);

    if (n > nparams)
	chisq /= (n - nparams);
    else
	chisq = 0;

    *p_chisq = chisq / (noise * noise);

    if (fit->method == SLOW_EXCHANGE_FIT)
	params_fit[1] = ABS(params_fit[1]);
    else if (fit->method == INVERSION_RECOVERY_FIT)
	params_fit[1] = params_fit[1] / BSCALE;

    if (params_dev)
    {
	for (i = 0; i < nparams; i++)
	    params_dev[i] = (float) sqrt((double) chisq*MAX(fit->covar[i][i], 0));
    }

    return  CCPN_OK;
}

void delete_fit(Fit_method fit)
{
    int method = fit->method;
    Log_linear_fit log_linear_fit;
    Non_linear_fit non_linear_fit;
    
    if (method == LOG_LINEAR_FIT)
    {
	log_linear_fit = (Log_linear_fit) fit;
        if (log_linear_fit->nalloc > 0)
            free_log_linear_memory(log_linear_fit);
	FREE(log_linear_fit, struct Log_linear_fit);
    }
    else if (method >= NONLINEAR2_FIT)
    {
	non_linear_fit = (Non_linear_fit) fit;
        if (non_linear_fit->nalloc > 0)
            free_nonlinear_memory(non_linear_fit);
	FREE(non_linear_fit, struct Non_linear_fit);
    }
    else
    {
	FREE(fit, struct Fit_method);
    }
}

CcpnStatus run_fit(Fit_method fit, int n, float *x, float *y, 
		float *params_fit, float *params_dev,
		float *y_fit, float *chisq, CcpnString error_msg)
{
    int method = fit->method;
    float noise = fit->noise;
    
    if (n < get_method_nparams(method))
    {
	COPY_VECTOR(y_fit, y, n);
	*chisq = 0;
	return  CCPN_OK;
    }

    if (method == LINEAR_FIT)
	return fit_linear(n, x, y, noise,
				params_fit, params_dev, y_fit, chisq, error_msg);
    else if (method == LOG_LINEAR_FIT)
	return fit_log_linear((Log_linear_fit) fit, n, x, y, noise,
				params_fit, params_dev, y_fit, chisq, error_msg);
    else if (method >= NONLINEAR2_FIT)
	return fit_nonlinear((Non_linear_fit) fit, n, x, y, noise,
				params_fit, params_dev, y_fit, chisq, error_msg);

    sprintf(error_msg, "unknown method %d", method);
    return CCPN_ERROR;
}

Fit_method new_fit(int method, float noise)
{
    Fit_method fit;
    Log_linear_fit log_linear_fit;
    Non_linear_fit non_linear_fit;

    if (method == LOG_LINEAR_FIT)
    {
	MALLOC_NEW(log_linear_fit, struct Log_linear_fit, 1);
	log_linear_fit->nalloc = 0;
	fit = (Fit_method) log_linear_fit;
    }
    else if (method >= NONLINEAR2_FIT)
    {
	MALLOC_NEW(non_linear_fit, struct Non_linear_fit, 1);
	non_linear_fit->nalloc = 0;
	fit = (Fit_method) non_linear_fit;
    }
    else /* ((method == NO_FIT) || (method == LINEAR_FIT)) */
    {
	MALLOC_NEW(fit, struct Fit_method, 1);
    }

    fit->method = method;
    fit->noise = noise;

    return fit;
}

static void sample_data(int n, float *x, float *x_dev, float *xw)
{
    int i;

    for (i = 0; i < n; i++)
	xw[i] = x[i] + x_dev[i]*normal01();
}

CcpnStatus run_fit_data(Fit_method fit, int niter, int n,
	float *x, float *y, float *x_dev, float *y_dev, float *xw, float *yw,
	float *params_fit, float *params_avg, float *params_dev,
	float *y_fit, float *chisq, CcpnString error_msg)
{
    int i, j, ngood, nparams, method = fit->method;
    double r;
    CcpnStatus status;

    nparams = get_method_nparams(method);
    ngood = 0;
    for (j = 0; j < nparams; j++)
	params_avg[j] = params_dev[j] = 0;

    for (i = 0; i < niter; i++)
    {
	sample_data(n, x, x_dev, xw);
	sample_data(n, y, y_dev, yw);
	status = run_fit(fit, n, xw, yw, params_fit, NULL, y_fit, chisq, error_msg);
	if (status == CCPN_OK)
	{
	    for (j = 0; j < nparams; j++)
	    {
		params_avg[j] += params_fit[j];
		params_dev[j] += params_fit[j] * params_fit[j];
	    }

	    ngood++;
	}
    }

    if (ngood < 2)
	RETURN_ERROR_MSG("not enough good fits when sampling");

    for (j = 0; j < nparams; j++)
    {
	params_avg[j] /= ngood;
	r = (double) (params_dev[j]/ngood - params_avg[j]*params_avg[j]);
	r = MAX(r, 0);
	params_dev[j] = (float) sqrt(r);
    }

    return run_fit(fit, n, x, y, params_fit, NULL, y_fit, chisq, error_msg);
}

static void bootstrap_data(int n, float *x, float *y, float *xw, float *yw)
{
    int i, k;

    for (i = 0; i < n; i++)
    {
	k = n * uniform01();
	k = MIN(k, n-1);
	xw[i] = x[k];
	yw[i] = y[k];
    }
}

/* Reference:
   Bootstrap Methods for Standard Errors, Confidence Intervals and Other
   Measures of Statistical Accuracy
   B. Efron and R. Tibshirani
   Statistical Science, 1986, Vol. 1, No. 1, 54-77
*/
CcpnStatus bootstrap_fit_data(Fit_method fit, int niter, int n,
	float *x, float *y, float *xw, float *yw,
	float *params_fit, float *params_avg, float *params_dev,
	float *y_fit, float *chisq, CcpnString error_msg)
{
    int i, j, ngood, nparams, method = fit->method;
    double r;
    CcpnStatus status;

    nparams = get_method_nparams(method);
    ngood = 0;
    for (j = 0; j < nparams; j++)
	params_avg[j] = params_dev[j] = 0;

    for (i = 0; i < niter; i++)
    {
	bootstrap_data(n, x, y, xw, yw);
	status = run_fit(fit, n, xw, yw, params_fit, NULL, y_fit, chisq, error_msg);
	if (status == CCPN_OK)
	{
	    for (j = 0; j < nparams; j++)
	    {
		params_avg[j] += params_fit[j];
		params_dev[j] += params_fit[j] * params_fit[j];
	    }

	    ngood++;
	}
    }

    if (ngood < 2)
	RETURN_ERROR_MSG("not enough good fits when sampling");

    for (j = 0; j < nparams; j++)
    {
	params_avg[j] /= ngood;
	r = (double) (params_dev[j] - ngood*params_avg[j]*params_avg[j]);
	r = MAX(r, 0);
	params_dev[j] = (float) sqrt(r / (ngood-1));
    }

    return run_fit(fit, n, x, y, params_fit, NULL, y_fit, chisq, error_msg);
}

int get_method_nparams(int method)
{
    static Bool printed = CCPN_FALSE;
    static int nparams[] = { 0, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 3, 3, 4, 4, 3, 2 };
    if (!printed && (sizeof(nparams)/sizeof(int) != NFIT_METHODS))
    {
        printf("INCONSISTENT nparams in get_method_params()\n");
        printed = CCPN_TRUE;
    }

    if ((method < 0) || (method >= NFIT_METHODS))
	return 0;

    return nparams[method];
}
