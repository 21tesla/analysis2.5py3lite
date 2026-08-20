/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*            Bayesian Inference / Massive Inference*/
/* */
/* Filename:  userstr.h*/
/* */
/* Purpose:   Define user structures.*/
/**/
/* History:   JS    25 Apr 2001 - 4 Feb 2003*/
/*            DJ O'D  5th July 2010                                            */
/*-----------------------------------------------------------------------------*/
/**/
#ifndef USERSTRH
#define USERSTRH

#include <Python.h>
#include "python_util.h"

/* CCPN for the strange block_file types in UserSpectrum */
#include "types.h"

typedef struct
{
    double  shape;
    int     pos_peaks;
    /* HEIGHT Cube[0] */
    double  max_height;
    double  min_height;
    double  dif_height;

    /* SIGMA Cube[1] */
    double  max_max_sigma;
    double  *max_sigma;
    double  *min_sigma;
    double  *dif_sigma;

    /* MEAN Cube[2] */
    double  *max_mean;
    double  *min_mean;
    double  *dif_mean;

    /* Q Cube vairies.. */
    double  max_Q;
    double  min_Q;
    double  dif_Q;
} UserPriorStr;

/* Data specific to the spectra files */
typedef struct
{
    char *  spc_file;
    int     ndim;
    int     real_ndim; /* for when we have SampledDataDim data */
    int     endianess;
    int  *  npoints;
    int     total_points;
    int  *  block_size; 

    int  *  sample_start;
    int  *  sample_end;
    int  *  sample_size;

/* block_file data - hopefully recoverable from the python world */
    int  *  dim_wrapped;
    int     bytes_per_point;        /* number of bytes per data point */
    Bool    big_endian;             /* whether data is big or little endian */
    Bool    padded;                 /* whether end blocks are padded */
    int     header;                 /* length of header in bytes */
    Bool    integer;                /* whether data is integer or real */
    Bool    writeable;              /* whether file is writeable or only readable */

    long double  MeanData;
} UserSpecStr;

typedef struct          /* COMMON PARAMETERS AND STATISTICS*/
{
    int       Ncell;      /* I   # object cells*/
    int       Nsample;    /*   O # output ensembles*/
    double    atoms;      /*   O <# atoms>*/
    double*   Mockbar;    /*   O <mock data>                        [Ndata]*/
    double*   PrPos;      /*   O Prob(object +ve)                   [Ncell]*/
    double*   Objbar;     /*   O Object average                     [Ncell]*/
    double*   Objdev;     /*   O Object std.dev.                    [Ncell]*/

    /* Custom additions! */
    double    Acc;
    float*    hist;

    int       try1;
    int       try2;
    int       insert1;
    int       insert2;
    int       delete1;
    int       empty;

    int       visual;

    int       Natoms;
    int       spec_dim;

    /* cont. work */
    double   *NormCube;
    double   *Mock;
    double   **model;
    int       MaxDimSize;
    void     *py_list;

    UserPriorStr *UserPrior;
    UserSpecStr  *UserSpec;

} UserCommonStr;

#endif /* USERSTRH */
