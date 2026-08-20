#include "fit1d.h"

#define EPS 1.0e-20
#define GOLD 1.618034
#define GLIMIT 100.0

#define TOL 0.03
#define C (GOLD-1)
#define R (1-C)

void bracket_minimum(float *p_ax, float *p_bx, float *p_cx, void *param, Fit1d_func func)
{
    float ax, bx, cx, fa, fb, fc, r, q, d, s, u, ulim, fu;

    ax = *p_ax;
    bx = *p_bx;

    fa = (*func)(ax, param);
    fb = (*func)(bx, param);

    if (fb > fa)
    {
	SWAP(ax, bx, float);
	SWAP(fa, fb, float);
    }

    cx = bx + GOLD*(bx-ax);
    fc = (*func)(cx, param);

    while (fb > fc)
    {
	r = (bx-ax) * (fb-fc);
	q = (bx-cx) * (fb-fa);
	d = q - r;
        s = d > 0? 1: -1;
	u = bx - 0.5*s*((bx-cx)*q-(bx-ax)*r)/MAX(ABS(d),EPS);
	ulim = bx + GLIMIT*(cx-bx);
	if ((bx-u)*(u-cx) > 0)
	{
	    fu = (*func)(u, param);
	    if (fu < fc)
	    {
		*p_ax = bx;
		*p_bx = u;
		*p_cx = cx;
	    }
	    else if (fu > fb)
	    {
		*p_ax = ax;
		*p_bx = bx;
		*p_cx = u;
	    }

	    u = cx + GOLD*(cx-bx);
	    fu = (*func)(u, param);
	}
	else if ((cx-u)*(u-ulim) > 0)
	{
	    fu = (*func)(u, param);
	    if (fu < fc)
	    {
		bx = cx;
		cx = u;
		u = cx + GOLD*(cx-bx);
	    }
	}
	else if ((u-ulim)*(ulim-cx) >= 0)
	{
	    u = ulim;
	    fu = (*func)(u, param);
	}
	else
	{
	    u = cx + GOLD*(cx-bx);
	    fu = (*func)(u, param);
	}

	ax = bx;
	bx = cx;
	cx = u;
	fa = fb;
	fb = fc;
	fc = fu;
    }

    *p_ax = ax;
    *p_bx = bx;
    *p_cx = cx;
}

float golden_search(float ax, float bx, float cx, void *param, Fit1d_func func, float *xmin)
{
    float x0, x1, x2, x3, f1, f2;

    x0 = ax;
    x3 = cx;
    if (ABS(cx-bx) > ABS(bx-ax))
    {
	x1 = bx;
	x2 = bx + C*(cx-bx);
    }
    else
    {
	x2 = bx;
	x1 = bx + C*(ax-bx);
    }

    f1 = (*func)(x1, param);
    f2 = (*func)(x2, param);

    while (ABS(x3-x0) > TOL*(ABS(x1)+ABS(x2)))
    {
	if (f2 < f1)
	{
	    x0 = x1;
	    x1 = x2;
	    x2 = R*x1 + C*x3;
	    f1 = f2;
	    f2 = (*func)(x2, param);
	}
	else
	{
	    x3 = x2;
	    x2 = x1;
	    x1 = R*x2 + C*x0;
	    f2 = f1;
	    f1 = (*func)(x1, param);
	}
    }

    if (f1 < f2)
    {
	*xmin = x1;
	return f1;
    }
    else
    {
	*xmin = x2;
	return f2;
    }
}

