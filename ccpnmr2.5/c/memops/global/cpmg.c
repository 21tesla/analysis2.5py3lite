#include "cpmg.h"

#include "fit1d.h"

/* for slow exchange, z is kex */
/* for fast exchange, z is kex - sqrt(kex^2 - dw^2) */
/* in both cases we are optimising over dw */
typedef struct _CpmgParam
{
    float R2max;
    float z;
    int n;
    float *x;
    float *y;
}  Cpmg_param;

#ifdef WIN64
double acosh(double x)
{
    if (x < 1)
        return 0.0; /* arbitrary */
      
    return log(x + sqrt(x*x-1));
}
#endif

static float diff2(int n, float *x, float *y, float *params)
{
    int i;
    float d, d2;

    d2 = 0;
    for (i = 0; i < n; i++)
    {
	d = cpmg3(x[i], params) - y[i];
	d2 += d * d;
    }

    return d2;
}

static float slow_exchange_func(float dw, void *param)
{
    Cpmg_param *p = (Cpmg_param *) param;
    int n = p->n;
    float params[3], *x = p->x, *y = p->y;

    params[0] = p->R2max;
    params[1] = p->z; /* kex */
    params[2] = dw;

    return diff2(n, x, y, params);
}

static float fast_exchange_func(float dw, void *param)
{
    Cpmg_param *p = (Cpmg_param *) param;
    int n = p->n;
    float params[3], *x = p->x, *y = p->y, z = p->z;

    params[0] = p->R2max;
    params[1] = 0.5 * (z*z + dw*dw) / z; /* kex */
    params[2] = dw;

    return diff2(n, x, y, params);
}

float cpmg3(float x, float *params)
{
/*
    see Mulder, Mittermaier, Hon, Dahlquist and Kay
    Nature Structural Biology 8 (2001) 932-935
    http://www.nature.com/nsmb/journal/v8/n11/pdf/nsb1101-932.pdf
    x = nu = 1/(2 tau)
*/
    float R2max = params[0], kex = params[1], dw = params[2];
    float psi, Dp, Dm, etap, etam, v;

    psi = kex*kex - dw*dw;
    if (psi > 0)
    {
        Dp = kex*kex / psi;
        Dm = dw*dw / psi;
        etap = sqrt(psi) / (2.0*x);
        v = Dp*cosh(etap) - Dm;
    }
    else
    {
        psi = -psi;
        Dp = dw*dw / psi;
        Dm = kex*kex / psi;
        etam = sqrt(psi) / (2.0*x);
        v = Dp - Dm*cos(etam);
    }

    return R2max + 0.5*kex - x*acosh(v);
}

void cpmg3_fast_init_params(int n, float *x, float *y, float *params_fit)
{
    int i;
    float ymin, ymax, z, R2max, ax, bx, cx, dw;
    Cpmg_param param;

    ymin = ymax = y[0];
    for (i = 1; i < n; i++)
    {
        ymax = MAX(ymax, y[i]);
        ymin = MIN(ymin, y[i]);
    }

    param.n = n;
    param.x = x;
    param.y = y;
    R2max = param.R2max = ymin;
    z = param.z = 2*(ymax-ymin);

    /* strangely, higher multiples of z mean lower dw/kex */
    ax = 8 * z;
    bx = 12 * z;
    bracket_minimum(&ax, &bx, &cx, &param, fast_exchange_func);

    golden_search(ax, bx, cx, &param, fast_exchange_func, &dw);

    params_fit[0] = R2max;
    params_fit[1] = 0.5 * (z*z + dw*dw) / z; /* kex */
    params_fit[2] = dw;
}

void cpmg3_slow_init_params(int n, float *x, float *y, float *params_fit)
{
    int i;
    float ymin, ymax, z, R2max, ax, bx, cx, dw;
    Cpmg_param param;

    ymin = ymax = y[0];
    for (i = 1; i < n; i++)
    {
        ymax = MAX(ymax, y[i]);
        ymin = MIN(ymin, y[i]);
    }

    param.n = n;
    param.x = x;
    param.y = y;
    R2max = param.R2max = ymin;
    z = param.z = 2*(ymax-ymin);

    ax = 4 * z;
    bx = 6 * z;
    bracket_minimum(&ax, &bx, &cx, &param, slow_exchange_func);

    golden_search(ax, bx, cx, &param, slow_exchange_func, &dw);

    params_fit[0] = R2max;
    params_fit[1] = z;
    params_fit[2] = dw;
}

float cpmg4(float x, float *params)
{
/*
    see Mulder, Mittermaier, Hon, Dahlquist and Kay
    Nature Structural Biology 8 (2001) 932-935
    http://www.nature.com/nsmb/journal/v8/n11/pdf/nsb1101-932.pdf
    x = nu = 1/(2 tau)
*/
    float R2max = params[0], kAB = params[1], kBA = params[2], dw = params[3];
    float psi, zeta, kex, Dp, Dm, etap, etam, v, t, s;

    kex = kAB + kBA;
    psi = kex*kex - dw*dw;
    zeta = 2*dw*(kAB - kBA);
    t = sqrt(psi*psi + zeta*zeta);
    s = (psi + 2*dw*dw) / t;

    Dp = 0.5 * (1 + s);
    Dm = 0.5 * (-1 + s);
    etap = sqrt(0.5*(psi+t)) / (2.0*x);
    etam = sqrt(0.5*(-psi+t)) / (2.0*x);
    v = Dp*cosh(etap) - Dm*cos(etam);

    return R2max + 0.5*kex - x*acosh(v);
}

static void cpmg4_init_params(float *params_fit, float *params3_fit)
{
    float pAB = 0.5, pBA = 1 - pAB;

    params_fit[0] = params3_fit[0]; /* R2max */
    params_fit[1] = pAB * params3_fit[1]; /* kAB = pAB * kex */
    params_fit[2] = pBA * params3_fit[1]; /* kBA = pBA * kex */
    params_fit[3] = params3_fit[2]; /* dw */
}

void cpmg4_fast_init_params(int n, float *x, float *y, float *params_fit)
{
    float params3_fit[3];
    cpmg3_fast_init_params(n, x, y, params3_fit);
    cpmg4_init_params(params_fit, params3_fit);
}

void cpmg4_slow_init_params(int n, float *x, float *y, float *params_fit)
{
    float params3_fit[3];
    cpmg3_slow_init_params(n, x, y, params3_fit);
    cpmg4_init_params(params_fit, params3_fit);
}

