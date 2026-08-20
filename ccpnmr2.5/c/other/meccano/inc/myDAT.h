#ifndef _MYDAT_H_
#define _MYDAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include "../inc/myMACRO.h"


typedef enum { CC_, CH_, NH_, CN_, C1H, C2H, CSA, CB_, HA_, PP_COUPL_NB } pepPlaCoupl_e;

typedef struct {
    int             couplNb;
    gsl_vector_int *res;
    gsl_vector     *y, *sigma;
} lgRgCoupl_t;                  /* Long range coupling */

typedef struct {
    gsl_vector     *y, *sigma;
    gsl_vector_int *flag;
    lgRgCoupl_t     HnHn;
    lgRgCoupl_t     HaHn;
    lgRgCoupl_t     CHn;
    lgRgCoupl_t     HnC;
    lgRgCoupl_t     CbHn;
    lgRgCoupl_t     HnCb;
} medium_t;

typedef struct {
    int             pepPlaNb;
    int             totalCouplNb;
    medium_t       *medium;
    double          csa11, csa22;
    int             csiFlag;
    double          phi0, psi0, phiEr, psiEr;
    int             phiFlag, psiFlag;
    lgRgCoupl_t     hbscNO, hbscON;
} pepPlaData_t;

typedef struct {
    int             firstPepPlaNb;
    int             totalPepPlaNb;
    int             mediaNb;
    int             hnHnFlag, haHnFlag, cHnFlag, hnCFlag, cbHnFlag, hnCbFlag, hbscFlag, phiPsiFlag;
    pepPlaData_t   *pepPla;
} data_t;

extern data_t  *InitData(int mediaNb, int totalPepPlaNb, int firstPepPlaNb);
extern void     DelData(data_t *data);

extern void     ReadPhiPsi(char *phiPsiFileName, data_t * data);
extern void     SetPhiPsi(int num, double r2, double r3, double r4, double r5, data_t *data);
extern void     ReadRdcFile(char *rdcFileName, int m, data_t * data);
extern void     SetRdcData(int m, int num1, char *atom1, int num2, char *atom2, double y, double sigma, data_t *data);
extern void     ReadCsiFile(char *csiFileName, data_t * data);

extern void     DelData(data_t * data);

#endif
