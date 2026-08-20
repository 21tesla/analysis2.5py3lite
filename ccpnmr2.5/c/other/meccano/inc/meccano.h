#ifndef _MECCANO_h_
#define _MECCANO_h_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>
#include "../inc/mySTRUCT.h"
#include "../inc/myPPGEO.h"
#include "../inc/myRDC.h"
#include "../inc/myDAT.h"
#include "../inc/myPARAMS.h"

extern void     RunMeccanoFromFiles(char *outputFileName, char *mediaFileName,
	    char *csiFileName, char *talosFileName, char *hbscFileName,
	    char *sequenceFileName, char *ramaDBFileName, int firstPPlaneFrag,
	    int lastPPlaneFrag, int ppNbMin, int minValueBest, int maxValueBest);

extern void     RunMeccanoFromParams(int firstPPlaneFrag, int lastPPlaneFrag,
	    int ppNbMin, int minValueBest, int maxValueBest, param_t *params);

extern param_t *NewParams();

extern void     DelParams(param_t *params);

#endif
