#ifndef _MYSTRUCT_h_
#define _MYSTRUCT_h_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>
#include "../inc/myGEOMETRY.h"
#include "../inc/myPPGEO.h"
#include "../inc/myDAT.h"

#define STRUCT_FORMAT "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s\n"
#define XYZ(vector) gsl_vector_get(vector, 0), gsl_vector_get(vector, 1), gsl_vector_get(vector, 2)

typedef enum { GLY, ALA, ARG, ASN, ASP, CYS, GLN, GLU, HIS, ILE, LEU, LYS, MET, PHE, SER, THR, TRP, TYR, VAL, PRO, PPR, PRP, NB_AC_AM_TYPE } acAmType_t;

typedef struct {
    int             resNb;
    gsl_vector     *n_, *ca, *c_, *h_, *o_, *cb, *ha, *xAxis, *yAxis, *zAxis;
    int             n_Flag, caFlag, c_Flag, h_Flag, o_Flag, cbFlag, haFlag;
    double          phi, psi, tet;
    double          first[3];
    acAmType_t      acAmType;
} acAm_t;

typedef struct {
    int             firstPepPlaNb;
    int             totalPepPlaNb;
    int             windowSize;
    int             windowOffset;
    double          tet0, tetEr, dTet0;
    int             tetFlag;
    acAm_t         *acAm;
} polyPeptide_t;

extern polyPeptide_t *InitPolyPeptide(int totalPepPlaNb, int firstPepPlaNb);

extern void     DelPolyPeptide(polyPeptide_t * polyP);

extern void     AddCbHa(int ppn, polyPeptide_t * polyP);

extern void     OrientatePepPlaFwd(const gsl_vector * x, int firstPepPlaFlag, int ppn, pPlane_t * pPlane, polyPeptide_t * polyP);

extern void     OrientatePepPlaRev(const gsl_vector * x, int firstPepPlaFlag, int ppn, pPlane_t * pPlane, polyPeptide_t * polyP);

extern void     ReadSequenceFile(char *fileName, polyPeptide_t * polyP);

extern void     SetSequence(polyPeptide_t * polyP, int seq_num, char *seq_name);

extern void     AcAmType2Name(acAmType_t acAmType, char name[4]);

extern void     AcAmName2Type(char name[4], acAmType_t * acAmType);

extern polyPeptide_t *ReadPdbFile(char *fileName);

extern void     WritePdbFile(char *fileName, polyPeptide_t * polyP);

#endif
