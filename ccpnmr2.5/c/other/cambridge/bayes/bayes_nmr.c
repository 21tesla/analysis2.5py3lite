/* =========================================== 
======================COPYRIGHT/LICENSE START==========================

bayes_nmr.c: Part of the Bayes program

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
/* = bayesNMR main file: generic N dimension = */
/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.
 * =========================================== */

#include "bayes_nmr.h"

#ifdef max
#undef max
#endif
#define max(a,b) a > b ? a : b

#define DIM_3_STRING    "%10.6f,\t%4d,\t%10.6f,\t%10.6f,\t%10.6f\n"
#define DIM_4_STRING    "%10.6f,\t%4d,\t%10.6f,\t%10.6f,\t%10.6f,\t%10.6f\n"
#define DIM_5_STRING    "%10.6f,\t%4d,\t%14.4f,\t%10.6f,\t%10.6f,\t%10.6f,\t%10.6f\n"
#define DIM_6_STRING    "%10.6f,\t%4d,\t%14.4f,\t%10.6f,\t%10.6f,\t%10.6f,\t%10.6f,%10.6f\n"

#ifdef SMALLLNUMBER
#undef SMALLLNUMBER
#endif
#define SMALLLNUMBER    1.0E-08

#ifdef MAX_DATA_RUN_POINTS
#undef MAX_DATA_RUN_POINTS
#endif
#define MAX_DATA_RUN_POINTS     65536  /* 2048 * 32 */

int bayesNMR(
    /* REGION provided by getRegion macro */
    char *  spc_file,
    int     ndim,
    int     endianess,
    int  *  npoints, 
    int  *  block_size, 
    int  *  sample_start,
    int  *  sample_end, 
    float   max_height,
    float   min_height,
    double  *   max_sigma,
    double  *   min_sigma,
    double  max_Q,
    double  min_Q,

    int  *  dim_wrapped,
    /* PARAMS user specified settings */
    int     shape,
    int     pos_peaks, 
    int     min_atoms,
    int     max_atoms,
    double  rate, 
    void  * py_list)
{
/* Addressing*/
    int             ENSEMBLE = 4;      /*  >= 1*/
    ObjectStr       *Objects;           /* [ENSEMBLE]*/
    CommonStr       Common[1];
    UserCommonStr   UserCommon[1];
    UserPriorStr    UserPrior[1];
    UserSpecStr     UserSpec[1];

    double          *Mockbar;           /* [Ndata]*/
    double          *PrPos;             /* [Ncell]*/
    double          *Objbar;            /* [Ncell]*/
    double          *Objdev;            /* [Ncell]*/
    double          *Accuracy;          /* [Ndata]*/

    int             i;
    int             * spc_npoints, * spc_nblock, * spc_dim_wrapped;
    int             * sample_size;

    int             Ncell, Ndata, real_ndim;
    float           * data_flt, data_max = 0.0;
    double          * data;
    double          * dif_sigma, * max_mean, * min_mean, * dif_mean;

    int             code, plane_flag = -1;
    double          max_max_sigma;
    double          Z;                  /* some stats value */

#ifdef __DEBUG__
/* varius timing stuff */
    time_t  t0, t1; /* time_t is defined on <time.h> and <sys/types.h> as long */
    clock_t c0, c1; /* clock_t is defined on <time.h> and <sys/types.h> as int */
#endif /*__DEBUG__*/

    /* If we're taking stuff on the plane, we're not really interested in the planes 
    orthogonal dims are we ? */
    real_ndim = ndim;
    if (sample_start[ndim-1] == sample_end[ndim-1])
    {
        plane_flag = sample_end[ndim-1];
        ndim--;
#ifdef __DEBUG__
        printf("*** Plane found : %d\n", plane_flag);
#endif /*__DEBUG__*/
    }

    /* alloc, set params and set into structures */
    spc_npoints = npoints;
    spc_nblock = block_size;
    spc_dim_wrapped = dim_wrapped;

    UserSpec->spc_file = spc_file;
    UserSpec->ndim = ndim;
    UserSpec->real_ndim = real_ndim;
    UserSpec->endianess = endianess;
    UserSpec->npoints = spc_npoints;
    UserSpec->block_size = spc_nblock;
    UserSpec->dim_wrapped = spc_dim_wrapped;


    MALLOC( sample_size, int, ndim);
    for( UserSpec->total_points = 1, i = 0; i < ndim; i++)
    {
        UserSpec->total_points *= spc_npoints[i];
        sample_size[i] = sample_end[i] - sample_start[i];

#ifdef __DEBUG__
        printf("spc_npoints %d\n", spc_npoints[i]);
        printf("*** Sample Size (dim %1d) : %6d\n", i, sample_size[i]);
        printf("*** UserSpec->total_points %d\n", UserSpec->total_points);
#endif /*__DEBUG__*/
    }

    UserSpec->sample_start = sample_start;
    UserSpec->sample_end   = sample_end;
    UserSpec->sample_size  = sample_size;

    for(Ndata = 1, i = 0; i < ndim; i++)
        Ndata *= sample_size[i];        /* Number of data points - ie size of x (or y) dim in spectrum*/
    Ncell = Ndata;                      /* Check... */

    if( Ndata > MAX_DATA_RUN_POINTS )
    {
        printf("&&& ERROR: Large spectral region chosen for Bayes separation (try something smaller)\n");
        printf("&&& Size chosen %d > limit %d\n", Ndata, MAX_DATA_RUN_POINTS);
        return 1;
    }

    /* for N-d + 1 spare (in case e_q)*/
    if ( ndim > 1)
    {
        MALLOC(UserCommon->NormCube, double, (2 * ndim) + 2);
    }
    else
    {
        MALLOC(UserCommon->NormCube, double, (2 * ndim) + 3);
    }
    MALLOC(UserCommon->Mock, double, Ndata);

    Common->MinAtoms = min_atoms;
    Common->MaxAtoms = max_atoms;
    UserPrior->shape = shape;
    UserPrior->pos_peaks = pos_peaks;

    /* set the priors - don't think this needs to be another function really */
    if (min_height == 0)
    {
        /*printf("*** min_height == 0 can cause problems, setting as 0.0001\n");*/
        min_height = 0.0001;
    }

    UserPrior->max_height =   max_height;
    UserPrior->min_height =   min_height;
    UserPrior->dif_height =   max_height - min_height;

    /* SIGMA Cube[1] */
    MALLOC(dif_sigma, double, ndim);

    for( max_max_sigma = 0, i = 0; i < ndim; i++ )
    {
        dif_sigma[i] = max_sigma[i] - min_sigma[i];
        max_max_sigma = max( max_max_sigma, max_sigma[i] );
    }

    UserPrior->max_sigma = max_sigma;
    UserPrior->min_sigma = min_sigma;
    UserPrior->dif_sigma = dif_sigma;
    UserPrior->max_max_sigma = max_max_sigma;

    /* MEAN Cube[2] */
    MALLOC(max_mean, double, ndim);
    MALLOC(min_mean, double, ndim);
    MALLOC(dif_mean, double, ndim);

    for( i = 0; i < ndim; i++ ) 
    {
        max_mean[i] = (double) sample_size[i];
        min_mean[i] = 0.0;
        dif_mean[i] = max_mean[i] - min_mean[i];
    }
    UserPrior->max_mean = max_mean;
    UserPrior->min_mean = min_mean;
    UserPrior->dif_mean = dif_mean;

    /* Q Cube vairies.. */
    UserPrior->max_Q = max_Q;
    UserPrior->min_Q = min_Q;
    UserPrior->dif_Q = max_Q - min_Q;

    /* get and set the data */
    if (get_sample_from_spec(UserSpec, &data_flt, plane_flag) != 0)
    {
        printf("&&& get_sample_from_spec: Error getting sample from spec, quitting\n");
        return 1;
    }
    if (UserSpec->total_points != Ndata)
    {
        printf("&&& UserSpec->total_points != Ndata: Error getting sample from spectrum\n");
        printf("&&& %d != %d\n", UserSpec->total_points, Ndata);
        return 1;
    }
    MALLOC(data, double, Ndata);
    for(i = 0; i < Ndata; i++)
    {
#ifdef WIN64
        if( _isnan( data_flt[i]) )
#else      
        if( isnan( data_flt[i]) )
#endif /*WIN32*/
        {
            printf("&&& Nan in input data\n");
            return 1;
        }
        
        data[i] = (double) data_flt[i];
        if (data_flt[i] > data_max) data_max = data_flt[i];
    }

#ifdef __DEBUG__
    printf("*** MAX DATA (%d) %f\n", Ndata, data_max);
#endif /*__DEBUG__*/

    FREE_TYPE(data_flt, float);

    /* set up one dim model shape arrays in UserCommon->model */
    {
        int max_dim_size = 0;
        for( i = 0; i < ndim; i++ )
            if ( max_dim_size < sample_size[i] )
                max_dim_size = sample_size[i];
        UserCommon->MaxDimSize = max_dim_size;
    }

    MALLOC2(UserCommon->model, double, ndim, UserCommon->MaxDimSize);

/* Initialise BayeSys*/
    MALLOC(Mockbar, double, Ndata);
    MALLOC(Objbar, double, Ncell);      /* Not used until after BayeSys so */
    MALLOC(Objdev, double, Ncell);      /*  something to do with the analysis*/
    MALLOC(PrPos, double, Ncell);

    MALLOC(Objects, ObjectStr, ENSEMBLE);

    Common->Alpha       = -1.0;         /* +ve for Poisson, -ve for geometric*/
    Common->ENSEMBLE    = ENSEMBLE;     /* # objects in ensemble*/
    Common->Method      =  0;           /* Algorithm method: all on -1, simplest 0, (no LS2 -3) - John recommended 0 for pp */
    Common->Rate        = rate;         /* Speed of calculation (dimensionless)*/
    Common->Iseed       = -4321;        /* Random seed, -ve is time seed*/

    MALLOC(Accuracy, double, Ndata);

    UserSpec->MeanData = 0;
    for(i = 0; i < Ndata; i++)
    {
        /* Initial guess - gets rescaled */
        Accuracy[i] = min_height;
        UserSpec->MeanData += (long double) data[i] /* * 1.0E-28*/;
    }
#ifdef WIN64
        if( _isnan(UserSpec->MeanData) )
#else      
        if( isnan(UserSpec->MeanData) )
#endif /*WIN32*/
    {
        printf("&&& MeanData is Nan :-|\n");
        return 1;
    }

    UserSpec->MeanData /= Ndata;

#ifdef __DEBUG__
    printf("\nTotal data = %f / Ndata %d\n", (double) UserSpec->MeanData, Ndata);
    printf("Mean data value = %f\n", (double) UserSpec->MeanData);
#endif /*__DEBUG__*/



/* Initialise Likelihood*/
    Common->Ndata       = Ndata;        /* Number of data points in real array*/
    Common->Data        = data;         /*  the real data*/
    Common->Acc         = Accuracy;     /* Accuracy set as noise threshold... ??? */

/* Initialise statistics*/
    Common->UserSpec    = (void*)UserSpec;
    Common->UserPrior   = (void*)UserPrior;
    Common->UserCommon  = (void*)UserCommon;	
    UserCommon->Ncell   = Ncell;      /* # cells*/
    UserCommon->Mockbar = Mockbar;    /* <mock data>                  [Ndata]*/
    UserCommon->PrPos   = PrPos;      /* Pr(cell occupied>            [Ncell]*/
    UserCommon->Objbar  = Objbar;     /* <intensity per cell>         [Ncell]*/
    UserCommon->Objdev  = Objdev;     /* std.dev. intensity per cell  [Ncell]*/
    UserCommon->Nsample = 0;          /* Nsamples from posterior */
    UserCommon->atoms   = 0.0;

/* Dan specific*/
    UserCommon->Acc     = 1.0;
    UserCommon->empty   = 0;
    UserCommon->try1    = 0;
    UserCommon->insert1 = 0;
    UserCommon->try2    = 0;
    UserCommon->insert2 = 0;
    UserCommon->delete1 = 0;

    UserCommon->spec_dim = ndim;

    if (shape == 5) Common->Ndim = (2 * ndim) + 2;      /* # of Cube dimensions (Unless using Wolfgang e_q)*/
    else            Common->Ndim = (2 * ndim) + 1;      /* # of Cube dimensions (Unless using Wolfgang e_q)*/

    Common->Valency     = 0;            /* # MassInf fluxes per atom*/

    for(i = 0; i < ENSEMBLE; i++ )      /* UserBuild uses mock data*/
        MALLOC(Objects[i].Mock, double, Ndata);

#ifdef __DEBUG__
    printf("*** Ndim: %d\n", Common->Ndim );
#endif /*__DEBUG__*/

/* open output files (pointers to) */

    UserCommon->py_list = py_list;

/* will be some more setting to be done, but moving on */
    for( i = 0; i < Ndata; i++ )
        Mockbar[i] = 0.0;
    for( i = 0; i < Ncell; i++ )
    {
        Objbar[i] = Objdev[i] = PrPos[i] = 0.0;
    }

#ifdef __DEBUG__
    printf("MaxAtoms:\t(>1)\t\t%d\nMinAtoms:\t(>1)\t\t%d\nENSEMBLE:\t(>1)\t\t%d\nRate:\t\t(>0)\t\t%5.3f\nNdim\t\t(>1)\t\t%d\n", \
            Common->MaxAtoms, Common->MinAtoms, Common->ENSEMBLE, Common->Rate, Common->Ndim);

    if (ndim > 1)
    {
        if ((UserPrior->dif_sigma[0] == 0) || (UserPrior->dif_sigma[1] == 0))
        {
            printf("bayesNMR: dif_sigma not correct: [0]: %f\t\t[1]: %f\n", UserPrior->dif_sigma[0], UserPrior->dif_sigma[1]);
            return 0;
        }
    }
    else if (UserPrior->dif_sigma[0] == 0)
    {
        printf("bayesNMR: dif_sigma not correct: [0]: %f\n", UserPrior->dif_sigma[0]);
        return 0;
    }
#endif /*__DEBUG__*/

/* Compute ...................................................................*/
    printf("\n*** Bayesian Peak Separator started ***\n");
    printf(  "***   N Peaks -   Cool      -    N Iter -     Likelihood\n");
#ifdef __DEBUG__
    t0 = time(NULL);
    c0 = clock();
#endif /*__DEBUG__*/
    code = BayeSys3(Common, Objects);
#ifdef __DEBUG__
    t1 = time(NULL);
    c1 = clock();
#endif /*__DEBUG__*/
    printf("\n***          Peak Separator ends ***\n");
/* ................................................................... compute*/
#ifdef __DEBUG__
    printf("Elapsed wall clock time:        %18ld\n", (long) (t1 - t0));
    printf("Elapsed CPU time:               %18.8f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
#endif /*__DEBUG__*/
/* some statistics of some variety */
    Z = Common->ENSEMBLE * UserCommon->Nsample; /* The same as out nsamp*/
    UserCommon->atoms /= Z;                                     /* <# atoms>*/
    for( i = 0; i < Ndata; i++ )
        Mockbar[i] /= Z;                                        /* <mock data>*/
    for( i = 0; i < Ncell; i++ )
    {
        PrPos[i] = 0.01*(double)(int)(0.5 + 100.0*PrPos[i]/Z);  /* Pr(positive)*/
        Objbar[i] /= Z;                                         /* <intensity>*/
        Objdev[i] /= Z;
        Objdev[i] -= Objbar[i] * Objbar[i];
        Objdev[i]  = (Objdev[i] > 0.0) ? sqrt(Objdev[i]) : 0.0; /* std.dev.*/
    }

/* Display results*/
#ifdef __DEBUG__
    printf("Random seed was   %d\n"  , Common->Iseed);
    printf("Success / CPU = %g / %g = %6.4f\n",
                    Common->Success, Common->CPU,
                    Common->Success / Common->CPU);
    printf("# ENSEMBLE samples was %d\n"  ,  UserCommon->Nsample);
    printf("< Atoms >   =%10.2f\n",          UserCommon->atoms);
    printf("Evidence    =%10.2f (log[e])\n", Common->Evidence);
    printf("Information =%10.2f (log[e])\n", Common->Information);
#endif /*__DEBUG__*/

    /* mimic python prompt */
    printf(">>> ");
    fflush(NULL);

    py_list = UserCommon->py_list;

    /* begin freeing up all the stuff */

    FREE_TYPE(UserCommon->NormCube, double);

    FREE_TYPE(sample_size, int);
    FREE_TYPE(dif_sigma, double);

    FREE_TYPE(max_mean, double);
    FREE_TYPE(min_mean, double);
    FREE_TYPE(dif_mean, double);

    FREE_TYPE(data, double);

    FREE_TYPE(Mockbar, double);
    FREE_TYPE(Objbar, double);
    FREE_TYPE(Objdev, double);
    FREE_TYPE(PrPos, double);

    FREE_TYPE(Objects, ObjectStr);
    FREE_TYPE(Accuracy, double);

    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Function:  UserMonitor*/
/**/
/* Purpose:   I have provided a new ensemble of objects.*/
/**/
/*        1.  For probabilistic exploration, restrict Common->cool to <= 1.*/
/*                    cool = 0 is prior, seen at the first call.*/
/*                0 < cool < 1 is "burn-in"*/
/*                    cool = 1 is evolve posterior,*/
/*                    cool -> infinity is maximum likelihood*/
/*            Return with a POSITIVE return code to exit BayeSys,*/
/*            usually after the end of an evolution stage with cool = 1.*/
/*            You can use Common->Nsystem, which counts calls to UserMonitor.*/
/**/
/*        2.  You should reset any nuisance parameters x in UserCommon and/or*/
/*            each of your UserObject structures, preferably by sampling from*/
/*            their conditonal probabilities Pr(x|Objects) and/or Pr(x|Object).*/
/**/
/*        3.  Collect your statistics and display your diagnostics. */
/*-----------------------------------------------------------------------------*/
/*  */
int UserMonitor(            /*   O  0 = continue, +ve = finish, -ve = abort*/
      CommonStr* Common,    /* I    Full general information*/
      ObjectStr* Objects)   /* I    Full ensemble of sample objects  [ENSEMBLE]*/
{
/* control*/
    int     Nrem;    /* # to go, as (2 * Nsample++) catches up with Nsystem++*/
/* desired statistics*/
    int     Ndata   = Common->Ndata;
    UserCommonStr* UserCommon = (UserCommonStr*)Common->UserCommon;

    int     Ncell   = UserCommon->Ncell;
    double* Objbar  = UserCommon->Objbar;  /* < intensity per cell >    [Ncell]*/
    double* Objdev  = UserCommon->Objdev;  /* std.dev.intensity per cell[Ncell]*/

/* local*/
    ObjectStr* Object;         /* individual object*/    
    UserPriorStr* prior = (UserPriorStr*) Common->UserPrior;

    double*    Cube;           /* individual atom coordinates in cube*/
    double     z;              /* gives individual atom flux*/
    double     flux;           /* cell-based flux from Cube[0] ##was Cube[1]*/
    int        cell;           /* cell number from Cube[2]     ##was Cube[0]*/
    int        k;              /* object counter*/
    int        r;              /* atom counter*/
    int        i,j;            /* general counter*/

    PyObject  * py_sample_obj;

    double      pr_unc, pr_llhood;
    int         pr_na;
    double      NormCube[MAX_NDIM];

    double      LhoodPrevious;

    Nrem = Common->Nsystem - (int) (2*UserCommon->Nsample);

/* Always print a line of diagnostics (cool, Nrem, logLikelihood[0])*/

    fprintf(stderr,                      /*(stderr flushes output immediately)*/

    "\t%2d    %8.6e \t %5d \t     %8.6e\r", \
          Objects[0].Natoms, Common->cool, Nrem, Objects[0].Lhood);


/*.............................................................................*/
/* Set any NUISANCE PARAMETERS you may have here*/
/* Could over-ride system's choices of FluxUnit (both Common and Objects) here*/
/* Could CALIBRATE data here (corrupting Evidence and Information)*/
/*.............................................................................*/

/* BURN-IN*/
    if( Common->cool < 1.0 )  /* final temperature not reached..*/
        return 0;             /* ...normal return, keep going*/
    else
        Common->cool = 1.0;

/*.............................................................................*/
/* EVOLUTION: Accumulate any desired statistics*/
    Object = &Objects[0];

    for( k = 0; k < Common->ENSEMBLE; k++ )
    {
        Object = &Objects[k];

        /* # atoms*/
        UserCommon->atoms += Object->Natoms;
        /* mock data*/

#ifdef _DONOTCOMPILE_
        for( j = 0; j < Ndata; j++ )
        {
            flux = 0.0;

            for( r = 0; r < Object->Natoms; r++ )
            {
                Cube = Object->Cubes[r];
                cell = (int)(Ncell * Cube[2]);
                if( cell == j )
                {
                    z = Cube[0];
                    flux += -10.0 * log(z);   /* direct BayeSys, q=10*/
                }
            }

            Objbar[j] += flux;                     /* sum intensity*/
            Objdev[j] += flux * flux;              /* sum intensity squared*/
        }
#endif /*_DONOTCOMPILE_*/

        for (j=0; j<Object->Natoms; j++)
        {
            /* use this here */

            cubes_priors_to_ncubes( Object->Cubes[j], prior, UserCommon->spec_dim, Common->Ndim, NormCube );

            /* set up print variables for each permutation */
            pr_unc    = Object->Cubes[j][Common->Ndim + Common->Valency];
            pr_llhood = Object->Lhood;
            pr_na     = Object->Natoms;

            /* NormCube structure:
                0   Height
                1   position    0
                2   sigma       0
               (3)  q
                3   position    1
                4   sigma       1
               (5)  q
                    etc.
            */

            /* print variables to file */

            py_sample_obj = (PyObject* ) PyList_New( Common->Ndim + 2 );

            PyList_SetItem(py_sample_obj, 0, PyFloat_FromDouble(    pr_llhood ) ); /* 0 */
            PyList_SetItem(py_sample_obj, 1, PyFloat_FromDouble( (int) pr_na  ) ); /* 1 */

            for( i = 0; i < 2*UserCommon->spec_dim + 1; i++ )
                PyList_SetItem(py_sample_obj, i+2, PyFloat_FromDouble( NormCube[i] ));

            if (Common->Ndim == ((UserCommon->spec_dim + 1) * 2) ) /* ie 4, 6, 8 etc*/
                PyList_SetItem(py_sample_obj, i+2, PyFloat_FromDouble( NormCube[UserCommon->spec_dim*2 + 1] ));

            if (PyList_Append( UserCommon->py_list, py_sample_obj) == -1)
            {
                printf("&&& PyList_Append fails\n");
            }
        } 
    }

    /* # samples of full ensemble*/
    UserCommon->Nsample++;

    if( Nrem > 0 )                                 /* not finished...*/
        return 0;                                  /* ..normal return*/

/*.............................................................................*/
/* EXIT WITH POSITIVE FINISH CODE WHEN EVOLUTION IS DEEMED COMPLETE*/
    return 1;
}

int get_sample_from_spec(UserSpecStr * UserSpec, float **data, int plane_flag)
{
    Block_file block_file;
    Mem_cache mem_cache;
    Line error_msg;
    struct  stat  statinfo;

    sprintf(error_msg, "getting spec sample from spectra");

    mem_cache = new_mem_cache(MEM_MAX_SIZE, NULL, NULL);

    if(stat(UserSpec->spc_file , &statinfo) == -1)
    {
        printf("&&& get_sample_from_spec: Error opening:\n %s\n", UserSpec->spc_file);
        return 1;
    }

    block_file = new_block_file(
        UserSpec->spc_file,     /* file */
        UserSpec->real_ndim,    /* ndim */
        UserSpec->npoints,      /* points */
        UserSpec->block_size,   /* block_size */
        UserSpec->dim_wrapped,  /* dim_wrapped ? */
        mem_cache,              /* mem_cache */
        BYTES_PER_WORD,         /* bytes_per_point (not 16 as in print out) */
        UserSpec->endianess,    /* big_endian */
        1,                      /* padded */
        0,                      /* header */
        0,                      /* integer */
        0,                      /* writeable */
        0);                     /* block_header */

    if (!block_file)
    {
        printf("&&& Error: Allocating block_file\n");
        return 1;
    }

    if (open_block_file(block_file) == CCPN_ERROR)
    {
        printf("&&& Error: opening block file '%s'", block_file->file);
        return 1;
    }

    {
        /* local hack for dodgy time data */
        int i, ndim = UserSpec->ndim;
        int temp_start[MAX_NDIM], temp_end[MAX_NDIM], temp_size[MAX_NDIM];
        int zeros[MAX_NDIM];

        for( i = 0; i < MAX_NDIM; i++ ) zeros[i] = 0;
#ifdef __DEBUG__
        printf("Total points %d : %d\n", UserSpec->total_points, plane_flag);
#endif /*__DEBUG__*/
        for( i = 0; i < ndim; i++ )
        {
            temp_start[i]   = UserSpec->sample_start[i];
            temp_end[i]     = UserSpec->sample_end[i];
            temp_size[i]    = temp_end[i] - temp_start[i];
#ifdef __DEBUG__
            printf("temp:       %4d      -> %4d\t", temp_start[i], temp_end[i]);
            printf("size:       %4d : plane %2d\n", temp_size[i], plane_flag);
#endif /*__DEBUG__*/
        }

        if (plane_flag > -1)
        {
            temp_start[ndim]   = plane_flag;
            temp_end[ndim]     = plane_flag + 1;
            temp_size[ndim]    = temp_end[ndim] - temp_start[ndim];
#ifdef __DEBUG__
            printf("temp:       %4d      -> %4d\t", temp_start[ndim], temp_end[ndim]);
            printf("size:       %4d : plane %2d\n", temp_size[ndim], plane_flag);
#endif /*__DEBUG__*/
        }

        if (get_box_block_file(block_file, &UserSpec->total_points, data, \
                            temp_start, temp_end, error_msg) == CCPN_ERROR)
        {
            printf("&&& Total points %d\n", UserSpec->total_points);
            printf("&&& %d %d\n", temp_start[0], temp_start[1]);
            printf("&&& %d %d\n", temp_end[0],  temp_end[1]);
            printf("&&& Error: get box from block file\n");
            return 1;
        }

#ifdef __DEBUG__
        for( i = 0; i < 5; i++ ) printf("%5d %f\n", i, (*data)[i] );
#endif /*__DEBUG__*/

#ifdef __DEBUG__
        printf("Total points %d\n", UserSpec->total_points);
        for( i = 0; i < ndim; i++ )
        {
            printf("%4d -> %4d\n", temp_start[i], temp_end[i]);
        }
        if (plane_flag > -1)
        {
            printf("%4d -> %4d\n", temp_start[ndim], temp_end[ndim]);
        }
#endif /*__DEBUG__*/
    }

    close_block_file(block_file);

    return 0;
}

void cubes_priors_to_ncubes( double *cube, UserPriorStr *prior, int ndim, int ComNdim, double *ncube )
{
    int i;

    /* Cube[0] : PEAK HEIGHT h */
    if (prior->pos_peaks == 1)       /* just search for positive peaks */
    {
        ncube[0] = prior->min_height   + cube[0] * prior->dif_height;
    }
    else                        /* include negative peaks */
    {
        if (cube[0] >= 0.5)
            ncube[0] =        (prior->min_height   + cube[0] * prior->dif_height);
        if (cube[0] < 0.5)
            ncube[0] =   -1 * (prior->min_height   + cube[0] * prior->dif_height);
    }

    /* cube[2*i+1] -> cube[2*i+1] : PEAK POSITION m */
    for ( i = 0; i < ndim; i++)
        ncube[2*i+1] = prior->min_mean[i]    + cube[2*i+1]   * prior->dif_mean[i];

    /* cube[2*i+2] -> cube[2*i+2] : SIGMA s */
    for ( i = 0; i < ndim; i++)
        ncube[2*i+2] = prior->min_sigma[i]   + cube[2*i+2]   * prior->dif_sigma[i];

    /* cube[3] : Q in 1d q - the Wolfgang extra variable */
    if( prior->shape == 5 )
        ncube[(ndim*2)+1]     = prior->min_Q          + cube[(ndim*2)+1]       * prior->dif_Q;

}

