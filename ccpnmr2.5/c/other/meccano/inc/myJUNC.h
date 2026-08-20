#ifndef _MYJUNC_H_
#define _MYJUNC_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../inc/myMACRO.h"
#include "../inc/myDAT.h"
#include "../inc/mySTRUCT.h"

extern double   CalcChi2_JuncGeo(polyPeptide_t * polyP, data_t * data, int resI);
extern void     DisplayJunctionParam(FILE * output, polyPeptide_t * polyP, int resI);
extern double   DisplayChi2_JuncGeo(FILE * output, polyPeptide_t * polyP, data_t * data, int resI);

#endif
