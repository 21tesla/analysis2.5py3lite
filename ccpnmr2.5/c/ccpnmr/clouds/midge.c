/*
======================COPYRIGHT/LICENSE START==========================

midge.c: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Alexei Grishaev, Miguel Llinas, Guillermo Bermejo, Wayne Boucher and Tim Stevens (Carnegie Mellon University and University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
M. Madrid, E. Llinas and M. Llinas (1991).
Model-Independent Refinement of Interproton Distances Generated from
H-1-NMR Overhauser Intensities. 
Journal Of Magnetic Resonance 93: 329-346.

A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/
#include "midge.h"

#include "diag_dbl.h"

/* convergence parameters */
#define  NUMERR  (1.0e-8)
#define  ALPHA  (1.0)
#define  CONVA  (1.0e-5)
#define  CONVR  (1.0e-5)
#define  DIAST  (1.0)
#define  EXPERR  (0.0)
#define  NEGEIG  (1.0e-5)

/* physical constants */
#define  HB    (1.05459e-27) /* hbar */
#define  GH    (2.6752e+4)   /* gamma for 1H */
#define  GN14  (1.9340e3)    /* gamma for 14N */
#define  GN15  (-2.7126e3)   /* gamma for 15N */
#define  GC    (6.7283e+3)   /* gamma for 13C */
#define  RNH   (0.96e-8)     /* H-N bond length */
#define  RCH   (1.070e-8)    /* H-C bond length */
#define  RGEM  (1.784e-8)    /* geminal bond length */
#define  TMET  (1.0e-10)     /* some time or other */
#define  CONST (GH*GH*GH*GH * HB*HB)

static CcpnStatus alloc_midge_memory(Midge midge)
{
    int n = midge->n;

    midge->copy = midge->aexp0 = midge->sp = NULL;
    midge->evalues = midge->work = NULL;

    MALLOC2(midge->copy, double, n, n);
    MALLOC2(midge->aexp0, double, n, n);
    MALLOC2(midge->sp, double, n, n);
    MALLOC(midge->evalues, double, n);
    MALLOC(midge->work, double, n);

    return CCPN_OK;
}

/*
static void print_vec(int n, double *vec, CcpnString msg)
{
    int i;

    printf("%s\n", msg);
    for (i = 0; i < n; i++)
	printf(" %8.3e", vec[i]);
    printf("\n");
}

static void print_vec2(int n, double *vec, int row)
{
    int i;

    printf("%d:", row);
    for (i = 0; i < n; i++)
	printf(" %8.3e", vec[i]);
    printf("\n");
}

static void print_mat(int n, double **mat, CcpnString msg)
{
    int i;

    printf("%s\n", msg);

    for (i = 0; i < n; i++)
	print_vec2(n, mat[i], i);
}
*/

static void olei_a(int n, double **amat, int *nhs)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	    amat[i][j] /= sqrt((double) (nhs[i] * nhs[j]));
    }
}

static void olei_r(int n, double **rmat, int *nhs)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	    rmat[i][j] *= sqrt((double) nhs[i]) / sqrt((double) nhs[j]);
    }
}

static void unolei_a(int n, double **amat, int *nhs)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	    amat[i][j] *= sqrt((double) (nhs[i] * nhs[j]));
    }
}

static void unolei_r(int n, double **rmat, int *nhs)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	    rmat[i][j] *= sqrt((double) nhs[j]) / sqrt((double) nhs[i]);
    }
}

static CcpnStatus diag(int n, double **matrix, double *evalues,
					double *work, CcpnString error_msg)
{
    return diagonalise_dbl(n, matrix, evalues, work, error_msg);
}

static void undiag(int n, double **output, double **input, double *eigenvalues)
{
    int i, j, k;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	{
	    output[i][j] = 0;

	    for (k = 0; k < n; k++)
		output[i][j] += input[i][k] * eigenvalues[k] * input[j][k];
	}
    }
}

static CcpnStatus a_to_r(int n, double *evalues, double tmix, CcpnString error_msg)
{
    int i, count = 0;

    for (i = 0; i < n; i++)
    {
	if (evalues[i] < 0)
	{
	    if (ABS(evalues[i]) < NEGEIG)
		evalues[i] = - evalues[i];
	    else
		evalues[i] = NEGEIG;

	    count++;
	}
    }

    if (count > 0)
    {
	if (count > 1)
	{
	    sprintf(error_msg, "%d negative eigenvalues\n", count);
	    return CCPN_ERROR;
	}
	else
	{
	    printf("%d negative eigenvalues\n", count);
	}
    }

    for (i = 0; i < n; i++)
	evalues[i] = - (double) log((double) evalues[i]) / tmix;

    return CCPN_OK;
}

static void r_to_a(int n, double *evalues, double tmix)
{
    int i;

    for (i = 0; i < n; i++)
	evalues[i] = (double) exp((double) (- tmix * evalues[i]));
}

static void cleanup_a(int n, double **amat)
{
    int i, j;
    double avg, minerr = 0; /* why 0? */

    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++) /* why i+1 and not 0? symmetric? */
	    minerr = MIN(minerr, amat[i][j]);
    }

/*
    printf("A numerical error = %8.3e\n", minerr);
*/

    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++)
	{
	    avg = 0.5 * (amat[i][j] + amat[j][i]);

	    if (avg < NUMERR)
		amat[i][j] = amat[j][i] = 0;
            else
		amat[i][j] = amat[j][i] = avg;
	}
    }
}

static void cleanup_r(int n, double **rmat)
{
    int i, j;
    double avg, maxerr = 0;

    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++) /* why i+1 and not 0? */
	    maxerr = MAX(maxerr, rmat[i][j]);
    }

/*
    printf("R numerical error = %8.3e\n", maxerr);
*/

    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++)
	{
	    avg = 0.5 * (rmat[i][j] + rmat[j][i]);

	    if (avg > NUMERR)
		rmat[i][j] = rmat[j][i] = 0;
            else
		rmat[i][j] = rmat[j][i] = avg;
	}
    }
}

static void calc_r_diag(int n, double **rmat, int *types,
	double rleak, double dipmet, double diparo, double dipch2,
	double dipnh, double dipch, double jself, double jcross, Bool c13_labelled)
{
    int i, j;
    double jrat = jself / jcross;

    for (i = 0; i < n; i++)
    {
	rmat[i][i] = rleak;

/* intramethyl-like relaxation */
	if ((types[i] == CH3_TYPE) || (types[i] == KHZ_TYPE))
	    rmat[i][i] += dipmet;

/* aromatic proton - pairs only */
	if (types[i] == HARO_TYPE)
	    rmat[i][i] += diparo;

/* no longer done
	if (types[i] == CH2_TYPE)
	    rmat[i][i] += dipch2;
*/

/* NH relaxation */
	if ((types[i] == HN_TYPE) || (types[i] == KHZ_TYPE))
	    rmat[i][i] += dipnh;

/* 13C relaxation (only in 13C-labelled sample) */
	if (c13_labelled && (types[i] != HN_TYPE))
	    rmat[i][i] += dipch;

	for (j = 0; j < n; j++)
	{
	    if ((i != j) && (rmat[i][j] != 0))
		rmat[i][i] += rmat[i][j] * jrat;
	}
    }
}

static Bool midgeconv(int n, double **aexp, double **sp, float *preverr)
{
    int i, j, err_cnt;
    double err;
    Bool conv;

    err = 0;
    err_cnt = 0;

    for (i = 0; i < n; i++)
    {
	for (j = 0; j < n; j++)
	{
	    if ((aexp[i][j] != 0) && (sp[i][j] != 0))
	    {
		err += ABS((aexp[i][j]-sp[i][j]) / aexp[i][j]);
		err_cnt++;
	    }
	}
    }

    err /= err_cnt;
    printf("error in A = %8.3e\n", err);
    printf("convergence in A = %8.3e\n", *preverr - err);

    if (ABS(err - *preverr) < CONVA)
	conv = CCPN_TRUE;
    else
	conv = CCPN_FALSE;

    *preverr = err;

/*
    printf("convergence = %s\n", conv ? "true" : "false");
*/

    return conv;
}

static void update_a(int n, double **amat, double **aexp, double **sp)
{
    int i, j;

/*
    diagonal and measured off-diagonal elements...
*/
    for (i = 0; i < n; i++)
    {
	for (j = i; j < n; j++) /* why is this i and not i+1? */
	{
	    amat[i][j] = aexp[i][j];
	    amat[j][i] = aexp[j][i];
	}

	amat[i][i] = ALPHA*sp[i][i] + (1.0-ALPHA)*amat[i][i];
    }

/*
    missing off-diagonal elements...
*/
    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++)
	{
	    if ((amat[i][j] == 0) && (sp[i][j] != 0))
	    {
		amat[i][j] = sp[i][j];
		amat[j][i] = sp[i][j];
	    }
	}
    }
}

Midge new_midge(int n, int *nhs, int *types)
{
    Midge midge;

    MALLOC_NEW(midge, struct Midge, 1);

    midge->n = n;
    midge->nhs = nhs;
    midge->types = types;

    if (alloc_midge_memory(midge) == CCPN_ERROR)
    {
	delete_midge(midge);
	return NULL;
    }

    return midge;
}

void delete_midge(Midge midge)
{
    int n = midge->n;

    FREE(midge->nhs, int);
    FREE(midge->types, int);
    FREE2(midge->copy, double, n);
    FREE2(midge->aexp0, double, n);
    FREE2(midge->sp, double, n);
    FREE(midge->evalues, double);
    FREE(midge->work, double);

    FREE(midge, struct Midge);
}

CcpnStatus run_midge(Midge midge, double **amat, double **rmat,
	int max_iter, float sf, float tmix, float tcor, float rleak,
	Bool n15_labelled, Bool c13_labelled,
	float *err, CcpnString error_msg)
{
    int i, j, cycle;
    Bool conv = CCPN_FALSE;
    int n = midge->n;
    int *nhs = midge->nhs;
    int *types = midge->types;
    double **copy = midge->copy;
    double **aexp0 = midge->aexp0;
    double **sp = midge->sp;
    double *evalues = midge->evalues;
    double *work = midge->work;

    double gn = n15_labelled ? GN15 : GN14;
    double tc = 1.0e-9 * tcor;
    double wh = 2.0 * PI * sf * 1.0e6; /* omega for H */
    double wn = wh * gn / GH; /* omega  for N */
    double wc = wh * GC / GH; /* omega for C */

/* correlation time (tc) - dependent params */
/* spectral density functions: j0,j1,j2 */
/* self relaxation rate constant: jself */
/* cross-realaxation rate constant: jcross */
    double j0 = CONST * tc;
    double j1 = CONST * tc / (1.0 + wh*wh*tc*tc);
    double j2 = CONST * tc / (1.0 + 4.0*wh*wh*tc*tc);
    double jself = 6.0*j2 + 3.0*j1 + j0;
    double jcross = 6.0*j2 - j0;

/* intramethyl relaxation rates */
/* Woessner, D., J. Chem. Phys. 36: 1 (1962) */
    double te = (tc * TMET) / (tc + TMET);
    double jmetself1 = (0.25*tc/(1.0+wh*wh*tc*tc) + (0.75*te/(1.0+wh*wh*te*te)))/pow(RGEM, 6.0);
    double jmetself2 = (0.25*tc/(1.0+4.0*wh*wh*tc*tc) + (0.75*te/(1.0+4.0*wh*wh*te*te)))/pow(RGEM, 6.0);
    double dipmet = 0.6 * CONST * (4.0*jmetself2 + jmetself1);

/* 14N-1H dipolar relaxation rates */
    double jn0 = tc / (1.0 + (wh-wn)*(wh-wn)*tc*tc);
    double jn1 = tc / (1.0 + wh*wh*tc*tc);
    double jn2 = tc / (1.0 + (wh+wn)*(wh+wn)*tc*tc);
    double dipnh = (1.0/10.0)*GH*GH*gn*gn*HB*HB *
	(jn0 + 3.0*jn1 + 6.0*jn2) / pow(RNH, 6.0);

/* 13C-1H dipolar relaxation rates */
    double jc0 = tc / (1.0 + (wh-wc)*(wh-wc)*tc*tc);
    double jc1 = tc / (1.0 + wh*wh*tc*tc);
    double jc2 = tc / (1.0 + (wh+wc)*(wh+wc)*tc*tc);
    double dipch = (1.0/10.0)*GH*GH*GC*GC*HB*HB *
	(jc0 + 3.0*jc1 + 6.0*jc2) / pow(RCH, 6.0);

/* aromatics and methylenes */
    double diparo = 0.3 * (4.0*j2 + j1) / pow(4.25e-8, 6.0); /* some bond length */
    double dipch2 = 0.3 * (4.0*j2 + j1) / pow(RGEM, 6.0);

    printf("***in midge***\n");
    printf("jself = %8.3e, jcross = %8.3e, jself/jcross = %8.3e\n",
				jself, jcross, jself/jcross);
    printf("te = %8.3e\n", te);
    printf("dipmet = %8.3e\n", dipmet);
    printf("dipnh = %8.3e\n", dipnh);
    printf("dipch = %8.3e\n", dipch);
    printf("diparo = %8.3e\n", diparo);
    printf("dipch2 = %8.3e\n", dipch2);

    tmix /= 1000;
    *err = 0;

    for (i = 0; i < n; i++)
    {
	for (j = i+1; j < n; j++)
	{
	    aexp0[i][j] = amat[i][j];
	    aexp0[j][i] = amat[j][i];
	}

	aexp0[i][i] = 0.0;
    }

/*
    print_mat(n, amat, "amat");
    print_mat(n, aexp0, "aexp0");
*/

    for (cycle = 0; (cycle < max_iter) && !conv; cycle++)
    {
	printf("midge cycle #%d\n", cycle);

	olei_a(n, amat, nhs);
/*
	print_mat(n, amat, "olei_a amat");
*/

	/* copy needed because diag overwrites matrix, and still need amat */
	for (i = 0; i < n; i++)
	{
	    for (j = 0; j < n; j++)
		copy[i][j] = amat[i][j];
	}

	sprintf(error_msg, "A diag: ");
	CHECK_STATUS(diag(n, copy, evalues, work, error_msg+strlen(error_msg)));
/*
	print_vec(n, evalues, "diag evalues");
	print_mat(n, copy, "diag copy");
*/

	sprintf(error_msg, "A to R: ");
	CHECK_STATUS(a_to_r(n, evalues, tmix, error_msg+strlen(error_msg)));
/*
	print_vec(n, evalues, "a_to_r evalues");
*/

	undiag(n, rmat, copy, evalues);
/*
	print_mat(n, rmat, "undiag rmat");
*/

	cleanup_r(n, rmat);
/*
	print_mat(n, rmat, "cleanup_r rmat");
*/

	unolei_r(n, rmat, nhs);
/*
	print_mat(n, rmat, "unolei_r rmat");
	print_mat(n, amat, "before unolei_a amat");
*/
	unolei_a(n, amat, nhs);
/*
	print_mat(n, amat, "unolei_a amat");
*/

	calc_r_diag(n, rmat, types, rleak, dipmet, diparo,
		dipch2, dipnh, dipch, jself, jcross, c13_labelled);
/*
	print_mat(n, rmat, "calc_r_diag rmat");
*/

	olei_r(n, rmat, nhs);
/*
	print_mat(n, rmat, "olei_r rmat");
*/
	/* copy needed because diag overwrites matrix, and still need rmat */
	for (i = 0; i < n; i++)
	{
	    for (j = 0; j < n; j++)
		copy[i][j] = rmat[i][j];
	}

	sprintf(error_msg, "R diag: ");
	CHECK_STATUS(diag(n, copy, evalues, work, error_msg+strlen(error_msg)));

	r_to_a(n, evalues, tmix);

	undiag(n, sp, copy, evalues);

	cleanup_a(n, sp);

	unolei_a(n, sp, nhs);

	conv = midgeconv(n, aexp0, sp, err);

	update_a(n, amat, aexp0, sp);
    }

    if (conv)
    {
	unolei_r(n, rmat, nhs);
    }
    else
    {
	printf(error_msg, "no convergence after %d iterations", max_iter);
	return CCPN_ERROR;
    }

    return CCPN_OK;
}
