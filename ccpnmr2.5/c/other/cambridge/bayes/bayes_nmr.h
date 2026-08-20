/* ============================== 
======================COPYRIGHT/LICENSE START==========================

bayes_nmr.h: Part of the Bayes program

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
/* = Header file for bayesNMR.c = */
/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.
 * ============================== */

#ifndef _incl_bayesNMR
#define _incl_bayesNMR

/* It seems both Marton Lowis (python) and C devs argue over 
   each others importance (C headers or Python headers first?)
*/
#include <Python.h>

#include <stdio.h>
#include <string.h>

/* timing stuff */
#include <time.h>
#include <sys/types.h>

#ifdef WIN64
#include <float.h>  /* for _isnan routine */
#endif

/* CCPN includes - mainly for get_sampl... */
#include <sys/stat.h>
#include "macros.h"
#include "block_file.h"
#include "mem_cache.h"

#include "bayesys3.h"
#include "userstr.h"

#include "userstr.h"

#ifndef MAX_NDIM
#define MAX_NDIM    5
#endif

#ifndef MAX_STRING
#define MAX_STRING      256
#endif


#define MEM_MAX_SIZE    134217728 /* 2 ^ 27 *//* poached from the ccpn running code */

#ifndef BYTES_PER_WORD
#define BYTES_PER_WORD 4
#endif

#ifdef  FREE_TYPE
#undef  FREE_TYPE
#endif
#define  FREE_TYPE(ptr, type) \
    {   if ((ptr) != (type *) NULL) \
        {   free((type *) (ptr));  (ptr) = (type *) NULL;   }   }


/* ARRAY ITERATION MACROS */
#define  INNER_PRODUCT(d, v1, v2, n) \
    {   int  I;  d = 0; \
        for (I = 0; I < (n); I++)  d += (v1)[I]*(v2)[I];   }

#define  INDEX_OF_ARRAY(index, array, cumul, n) \
    {   INNER_PRODUCT(index, array, cumul, n);  }

#define  ARRAY_OF_INDEX(array, index, cumul, n) \
    {   int I; long Ind = index; \
        for (I = (n)-1; I >= 0; I--) \
        {   array[I] = Ind / cumul[I];  Ind %= cumul[I];   }   }

#define  CUMULATIVE(cumul, array, total, n) \
    {   int I;  total = 1; \
        for (I = 0; I < n; I++) \
        {   (cumul)[I] = total;  total *= (array)[I];   }   }


int get_sample_from_spec(UserSpecStr * UserSpec, float ** data, int plane_flag);

/* given cubes, priors and shape return norm cubes */
void cubes_priors_to_ncubes( double *cube, UserPriorStr *prior, int ndim, int ComNdim, double *ncube );

extern int bayesNMR(
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
    void  * py_list);

/* Not really necessary - for the python wrappers */

typedef struct Bayes
{
     /* REGION provided by getRegion macro */
    char *  spc_file;
    int     ndim;
    int  *  npoints; 
    int  *  block_size; 
    int  *  sample_start;
    int  *  sample_end;
    int  *  sample_size;
    float   max_height;
    float   min_height;
    double  *   max_sigma;
    double  *   min_sigma;
    double  max_Q;
    double  min_Q;

    int  *  dim_wrapped;
    /* PARAMS user specified settings */
    int     shape;
    int     pos_peaks;
    int     min_atoms;
    int     max_atoms;
    double  rate;
}   *Bayes;

typedef struct {
  float    x;
  float    y;
  float    xerr;
  float    yerr;
} point;

#endif /* _incl_bayesNMR */

