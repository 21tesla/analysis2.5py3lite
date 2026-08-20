#ifndef _MYPARAMS_H_
#define _MYPARAMS_H_

#include "../inc/mySTRUCT.h"
#include "../inc/myRDC.h"
#include "../inc/myDAT.h"
#include "../inc/myPPGEO.h"
#include "../inc/myRAMACHANDRAN.h"

typedef struct {
    int             ppi;
    int             pPlaneNbInMin;
    int             windowStart;
    int             iterNb;
    pPlane_t       *pPlaneRef;
    tensor_t       *tensorStat;
    tensor_t       *tensorDyn;
    data_t         *data;
    polyPeptide_t  *polyP;
    ramaDB_t       *ramaDB;
    double          (*Func) (const gsl_vector *, void *);
} param_t;

#endif
