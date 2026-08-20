#ifndef _MYRDC_h_
#define _MYRDC_h_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "../inc/myMACRO.h"
#include "../inc/myGEOMETRY.h"
#include "../inc/myPPGEO.h"
#include "../inc/mySTRUCT.h"
#include "../inc/myDAT.h"

#define GH      267.669e6
#define GC      67.3027e6
#define GN      -27.125e6

typedef struct {
    double          Aa;
    double          Ar;
    double          alpha;
    double          beta;
    double          gamma;
    gsl_vector     *xCol;
    gsl_vector     *zCol;
    gsl_matrix     *rotMat;
    char            name[100];
} tensor_t;

extern tensor_t *InitTensor(int mediaNb);
extern void     ReadTensorFile(char *tensorName, tensor_t * tensor);
extern void     SetTensor(char *tensorName, double Aa, double Ar, double alpha, double beta, double gamma, tensor_t * tensor);
extern void     PrintTensor(FILE * output, tensor_t * tensor);
extern void     TotalPepPlaNbDet(char **media, int mediaNb, int *totalPepPlaNb, int *firstPepPlaNb);
extern void     UpdatePepPlaNbDet(int num1, char *atom1, int num2, char *atom2, int *fPepPla, int *lPepPla);
extern double   RdcStat(gsl_vector * vec, tensor_t * tensor);
extern double   RdcDynS(double S, const gsl_vector * vec, const tensor_t * tensor);
extern double   RdcDynGaf(double sig, const gsl_vector * axis, const gsl_vector * vec, const tensor_t * tensor);
extern double   Gaf1D(double Aa, double Ar, double ta, double pa, double t, double p, double sig);

extern double   CalcChi2Stat_PPRdc(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_PPRdc_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_CbHaRdc(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_CbHaRdc_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HHRdc_Rev(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HaHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HaHRdc_Rev(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HHRdc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HHRdc_Rev_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HaHRdc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_HaHRdc_Rev_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_CbHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   CalcChi2Stat_CHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);

extern double   DisplayChi2_PPRdc(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_CbHaRdc(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_HHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_HHRdc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_HaHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_HaHRdc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_CHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);
extern double   DisplayChi2_CbHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn);

extern void     DelTensor(tensor_t *tensor, int mediaNb);

#endif
