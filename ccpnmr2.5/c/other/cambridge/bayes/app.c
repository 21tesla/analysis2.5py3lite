/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
======================COPYRIGHT/LICENSE START==========================

app.c: Part of the Bayes program

Copyright (C) 2003-2010 Daniel O'Donovan (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the author: djo35@bioc.cam.ac.uk
=======================================================================
===========================REFERENCE START=============================
===========================REFERENCE END===============================

*/
/*            Bayesian Inference                                              */
/*                                                                            */
/* Filename:  bayesapp.c                                                      */
/*                                                                            */
/* Purpose:   Test operation of BayeSys with                                  */
/*            UserEmpty,UserTry1,UserTry2,UserInsert1,UserInsert2,UserDelete1.*/
/*============================================================================*/
/*                                                                            */
/* BayeSys3 can operate with arbitrary likelihood functions.                  */
/* In this example, an atom has 2 coordinates interpreted as                  */
/*            location in [0,1] = Cube[0]  (restricted to [0.0.75])           */
/* and                                                                        */
/*         flux in [0,infinity) = -10.0*log(Cube[1])                          */
/* to correspond with testing the operation of MassInf with FluxUnit0=10      */
/*                                                                            */
/* These are only example programs, involving inefficient calculations that   */
/* repeatedly build an object from scratch with "UserBuild" instructions.     */
/*                                                                            */
/*============================================================================*/
/* History:   JS       2 Jan 2002, 4 Jan 2003, 10 Feb 2003, 20 Aug 2003       */
/*            DJ O'D   2008-09-22.                                            */
/*            Copyright (c) 2002,2003 Maximum Entropy Data Consultants Ltd.   */
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bayesys3.h"
#include "userstr.h"
#include "bayes_nmr.h" /* only for MAX_NDIM */
#include "app.h"
#include "distribution.h"

/* CCPN */
#include "macros.h"

/* local functions */
static void generate_model(double * data, double ** model, int ndim, \
                            int * size, int ndata, \
                            double * param, int shape, double max_max_sigma);

static void generate_spectrum( double * data, double ** model, int ndim, \
                                    int * size, int ndata, \
                                    double * param, int shape, double max_max_sigma);


static int in_region( int ndim, double * param, int * position, double delta );

/*=============================================================================*/
/*                           SIMPLE CODE                                       */
/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserEmpty                                                        */
/*                                                                             */
/* Purpose:   Set Lhood = logLikelihood(coords) for empty sample               */
/*            with 0 atoms in Object, and initialise its other information.    */
/*                                                                             */
/*            Lhood := log L(empty)                                            */
/*                                                                             */
/*            I have already put the number 0 in Object->Natoms,               */
/*            so you can use a simple "UserBuild" instruction.                 */
/*-----------------------------------------------------------------------------*/
/* */
int UserEmpty(        /*   O  >=0, or -ve error code*/
double*    Lhood,     /*   O  loglikelihood*/
CommonStr* Common,    /* I O  general information*/
ObjectStr* Object)    /* I O  sample object, output new Lhood*/
{
    UserCommonStr *UserCommon = Common->UserCommon;

    UserBuild(Lhood, Common, Object, 0);

    UserCommon->empty++;
    return 1;
}

/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserTry1                                                         */
/*                                                                             */
/* Purpose:   d(logLikelihood(coords)) after supposedly adding one new atom    */
/*            to Object, without any side effects on it.                       */
/*                                                                             */
/*            If Valency = 0,                                                  */
/*              dLtry := d logL(...,x)                                         */
/*            else mimic the above with                                cool    */
/*              dLtry := (1/cool) log INTEGRAL dPrior(z) dlogL(...,x,z)        */
/*                                                                             */
/*            I have already put the new atom x after the last location        */
/*            of the Atoms list, at Object->Cubes[ Object->Natoms ], so        */
/*            if Valency=0 you can use simple "UserBuild" instructions.        */
/*-----------------------------------------------------------------------------*/
/* */
int UserTry1(         /*   O  +ve = OK, 0 = DO NOT USE, -ve = error*/
double*    dLtry,     /*   O  trial d(logLikelihood) value*/
CommonStr* Common,    /* I    general information*/
ObjectStr* Object)    /* I    sample object (DO NOT UPDATE Lhood)*/
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                      /* existing Lhood*/
    double Ltry;
    int    OK;

    UserCommonStr *UserCommon = Common->UserCommon;

    *dLtry = 0.0;

    OK = UserBuild(&Ltry, Common, Object, Natoms+1); /* trial*/


    if( OK > 0 )
        *dLtry = Ltry - Lhood;                          /* increment*/
    UserCommon->try1++;
    return  OK;                                         /* OK?*/
}

/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserTry2                                                         */
/*                                                                             */
/* Purpose:   d(logLikelihood(coords)) after supposedly adding one,            */
/*            then two new atoms to Object, without any side effects on it.    */
/*                                                                             */
/*   If Valency = 0,                                                           */
/*     dLtry1 := d logL(...,x1)                                                */
/*                                                                             */
/*     dLtry2 := d logL(...,x1,x2)                                             */
/*                                                                             */
/*   else mimic the above with                                    cool         */
/*     dLtry1 := (1/cool) log INTEGRAL dPrior(z1) dlogL(...,x1,z1)             */
/*                                                                         cool*/
/*     dLtry2 := (1/cool) log INTEGRAL dPrior(z1,z2) dlogL(...,x1,z1,x2,z2)    */
/*                                                                             */
/*            I have already put the new atoms x1,x2 after the last location   */
/*            of the Atoms list, at Object->Cubes[ Object->Natoms ] and        */
/*                               at Object->Cubes[ Object->Natoms + 1 ]        */
/*            so if Valency=0 you can use simple "UserBuild" instructions.     */
/*-----------------------------------------------------------------------------*/

int UserTry2(         /*   O  +ve = OK, 0 = DO NOT USE, -ve = error*/
double*    dLtry1,    /*   O  trial d(logLikelihood) value for 1st atom*/
double*    dLtry2,    /*   O  trial d(logLikelihood) value for both atoms*/
CommonStr* Common,    /* I    general information*/
ObjectStr* Object)    /* I    sample object (DO NOT UPDATE Lhood)*/
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                              /* existing Lhood*/
    double Ltry1, Ltry2;
    int    OK;

    UserCommonStr *UserCommon = Common->UserCommon;

    *dLtry1 = *dLtry2 = 0.0;

    OK = UserBuild(&Ltry1, Common, Object, Natoms+1);           /*trial for 1 more*/

    if( OK > 0 )
    {
        *dLtry1 = Ltry1 - Lhood;                                /* increment for 1*/

        OK = UserBuild(&Ltry2, Common, Object, Natoms+2);       /*trial for 2 more*/

        if( OK > 0 )
            *dLtry2 = Ltry2 - Lhood;                            /* increment for 2*/
    }
    UserCommon->try2++;
    return  OK;                             /* return OK only if both trials OK*/
}

/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserInsert1                                                      */
/*                                                                             */
/* Purpose:   Insert 1 new atom into Object, keeping it up-to-date, and        */
/*            set d(loglikelihood(coords)).                                    */
/*                                                                             */
/*            If Valency = 0,                                                  */
/*              dL := d logL(...,x)                                            */
/*                                                                             */
/*            else                                                             */
/*                                                cool                         */
/*              sample z from  Prior(z) L(...,x,z)                             */
/*                                                                             */
/*              and set    dL := d logL(...,x,z)  at the sampled z.            */
/*                                                                             */
/*            I have already put the new atom x at the last location of        */
/*            the updated Atoms list, at Object->Cubes[ Object->Natoms - 1 ],  */
/*-----------------------------------------------------------------------------*/

int UserInsert1(      /*   O  >=0, or -ve error code*/
double*    dL,        /*   O  d(loglikelihood)*/
CommonStr* Common,    /* I    general information*/
ObjectStr* Object)    /* I O  sample object*/
{
    int    Natoms = Object->Natoms;              /* new number*/
    double Lold   = Object->Lhood;               /* not yet updated*/
    double Lnew;

    UserCommonStr *UserCommon = Common->UserCommon;

    UserBuild(&Lnew, Common, Object, Natoms);       /* new updated state*/

    *dL = Lnew - Lold;                              /* I will update Object->Lhood*/
    UserCommon->insert1++;
    return 1;
}

/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserInsert2                                                      */
/*                                                                             */
/* Purpose:   Insert 2 new atoms into Object, keeping it up-to-date, and       */
/*            set d(loglikelihood(fluxes)).                                    */
/*                                                                             */
/*            If Valency = 0,                                                  */
/*              dL := d logL(...,x1,x2)                                        */
/*                                                                             */
/*            else                                                             */
/*                                                                cool         */
/*              sample z1,z2 from  Prior(z1,z2) L(...,x1,z1,x2,z2)             */
/*                                                                             */
/*              and set  dL := d logL(...,x1,z1,x2,z2)  at the sampled z1,z2.  */
/*                                                                             */
/*            I have already put the new atoms x1,x2 at the last location      */
/*            of the Atoms list, at Object->Cubes[ Object->Natoms - 2 ] and    */
/*                               at Object->Cubes[ Object->Natoms - 1 ]        */
/*            so if Valency=0 you can use a simple "UserBuild" instruction.    */
/*-----------------------------------------------------------------------------*/

int UserInsert2(      /*   O  >=0, or -ve error code*/
double*    dL,        /*   O  d(loglikelihood)*/
CommonStr* Common,    /* I    general information*/
ObjectStr* Object)    /* I O  sample object*/
{
    int    Natoms = Object->Natoms;              /* new number*/
    double Lold   = Object->Lhood;               /* not yet updated*/
    double Lnew;

    UserCommonStr *UserCommon = Common->UserCommon;

    UserBuild(&Lnew, Common, Object, Natoms);    /* new updated state*/

    *dL = Lnew - Lold;                           /* I will update Object->Lhood*/
    UserCommon->insert2++;
    return 1;
}

/* -= Not changed much from Steve's Code =-                                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserDelete1                                                      */
/*                                                                             */
/* Purpose:   Delete 1 old atom from Object, keeping it up-to-date, and        */
/*            set d(loglikelihood(fluxes)).                                    */
/*                                                                             */
/*            dL := d logL(...)                                                */
/*                                                                             */
/*            I have already put the old atom after the last location of       */
/*            the updated Atoms list, at Object->Cubes[ Object->Natoms ],      */
/*            so if Valency=0 you can use a simple "UserBuild" instruction.    */
/*-----------------------------------------------------------------------------*/

int UserDelete1(      /*   O  >=0, or -ve error code*/
double*    dL,        /*   O  d(loglikelihood)*/
CommonStr* Common,    /* I    general information*/
ObjectStr* Object)    /* I O  sample object*/
{
    int    Natoms = Object->Natoms;              /* new number*/
    double Lold   = Object->Lhood;               /* not yet updated*/
    double Lnew;

    UserCommonStr *UserCommon = Common->UserCommon;

    UserBuild(&Lnew, Common, Object, Natoms);    /* new updated state*/

    *dL = Lnew - Lold;                           /* I will update Object->Lhood*/
    UserCommon->delete1++;
    return 1;
}

/*=============================================================================*/
/*       Dummy procedure to link to BayeSys3 when not running MassInf          */
/*=============================================================================*/
int UserFoot(
double*    Cube,
CommonStr* Common,
int*       ibits,
double*    zbits,
int*       nbits)
{
    return (Cube && Common && ibits && zbits && nbits) ? 1 : 1;
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserBuild                                                        */
/*                                                                             */
/* Purpose:   Build Lhood, and any other data info such as Mock from scratch.  */
/*            Can be called with any number of atoms                           */
/*                                                                             */
/*            This example implementation has                                  */
/*            Mock[0] = SUM[atoms] flux ,             flux = - 10*log(Cube[1]) */
/*            Mock[1] = SUM[atoms] flux * location ,  location = Cube[0]       */
/*-----------------------------------------------------------------------------*/

int UserBuild(        /*   O  +ve = OK, 0 = DO NOT USE, -ve = error*/
double*    Lhood,     /*   O  loglikelihood */
CommonStr* Common,    /* I    General information*/
ObjectStr* Object,    /* I(O) Cubes in, perhaps Mock out*/
int        Natoms)    /* I    # atoms */
{
    double** Cubes  = Object->Cubes;  /* I    Cubes in [0,1)  [Natoms][Ndim]*/
    int      Ndata  = Common->Ndata;  /* I    # data*/
    double*  Data   = Common->Data;   /* I    Data            [Ndata]*/

    double*  Acc    = Common->Acc;   /* I    Data            [Ndata]*/

    UserCommonStr *UserCommon   = (UserCommonStr *) Common->UserCommon;
    UserPriorStr *prior         = (UserPriorStr *) Common->UserPrior; 
    UserSpecStr *spec           = (UserSpecStr *) Common->UserSpec;

    double      * Cube, * NormCube, * Mock, ** Model;
    int         i, k, ndim;
    double      mSum = 0, Lnorm = 0;
    const double Sqrt2Pi = 2.50662827463100050240;

    ndim = spec->ndim;
    Model = UserCommon->model;

    if (Natoms == 0)
    {
        *Lhood = 0;
        return 1;
    }

    /* make sure this is a zero array (ita better way?) */
    Mock = UserCommon->Mock;
    for(i = 0; i < Ndata; i++)
        Mock[i] = 0.0;

    NormCube = UserCommon->NormCube;
    for(i = 0; i < Common->Ndim; i++)
        NormCube[i] = 0.0;

    for(i = 0; i < ndim; i++)
    {
        for(k = 0; k < UserCommon->MaxDimSize; k++)
        {
            Model[i][k] = 0.0;
        }
    }

    
    /* Accumulate logLikelihood*/
    for( k = 0; k < Natoms; k++ )
    {
        Cube = Cubes[k];

        cubes_priors_to_ncubes( Cube, prior, ndim, Common->Ndim, NormCube );

        /* now get the data */
        generate_model(Mock, Model, spec->ndim, spec->sample_size, Ndata, NormCube, \
                            (int) prior->shape, prior->max_max_sigma);
    }

    generate_spectrum(Mock, Model, spec->ndim, spec->sample_size, Ndata, NormCube, \
                         (int) prior->shape, prior->max_max_sigma);

    mSum = 0;
    #pragma omp parallel \
        default( none ) \
        private( k ) \
        shared(  Ndata, Mock, Data, Acc ) \
        reduction(+:mSum, Lnorm)
    #pragma omp for
    for(k = 0; k < Ndata; k++)
    {
        mSum  += (Mock[k] - Data[k]) * Acc[k] * Acc[k] * (Mock[k] - Data[k]);
        Lnorm += log(Acc[k] / Sqrt2Pi);     // (fussy to include this)
    }

    *Lhood = Lnorm - mSum / Ndata;

    for( k = 0; k < Ndata; k++)
        if (Acc[k] > 0)
            Acc[k] = Acc[k] * sqrt( (float)Ndata / mSum );


    Object->Mock = Mock;
    UserCommon->Mockbar = Mock;
    UserCommon->Mock = Mock;

    return 1;
}

void generate_model( double * data, double ** model, int ndim, int * size, int ndata, \
                            double * param, int shape, double max_max_sigma)
{
    int i, j, n;

    double h;

    /* create one d arrays for each model */
    for( n = 0; n < ndim; n++ )
    {
        if (n == 0)
            h = param[0];
        else
            h = 1.;

        for( i = 0; i < size[n]; i++ )
        {
            /* Gaussian(double h, double s, double m, double x) */
            if (shape == 3) /* 3 : 2d Gauss */
                model[n][i] += Gaussian(     h, param[n*2+2], param[n*2+1], i);
            if (shape == 4) /* 4 : 2d Lorentz */
                model[n][i] += Lorentzian(   h, param[n*2+2], param[n*2+1], i);
            if (shape == 5) /* 4 : 2d Wolfgang */
                model[n][i] += Wolfgang(     h, param[n*2+2], param[n*2+1], param[n*2+3], i);
        }
    }

    return;
}

void generate_spectrum( double * data, double ** model, int ndim, int * size, int ndata, \
                                    double * param, int shape, double max_max_sigma)
{
    int i, j, n;
    int cumul_ar[MAX_NDIM], array[MAX_NDIM];
    double new_datum;


    CUMULATIVE(cumul_ar, size, ndata, ndim);

    for( i = 0; i < ndata; i++ )
    {

        ARRAY_OF_INDEX( array, i, cumul_ar, ndim);
        new_datum = 1.;

        for( n = 0; n < ndim; n++ )
        {
            new_datum *= model[n][array[n]];
        }

        data[i] += new_datum;
    }

    return;
}

/* If point is close to mean, calculate return 1, else return 0 */
int in_region( int ndim, double * param, int * position, double delta )
{
    int i, good = 1;
    double DELTA = delta * 2;

    for( i = 0; i < ndim; i++ )
    {
        /* if position fails set as bad */
        if (! (( param[i+1] - DELTA * param[i+2] <= position[i] ) && ( position[i] <= param[i+1] + DELTA * param[i+2] )) )
        {
            good = 0;
            break;
        }
    }

    return good;
}


