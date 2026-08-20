#ifndef _MYHBSC_H_
#define _MYHBSC_H_

#include "../inc/myMACRO.h"
#include "../inc/mySTRUCT.h"
#include "../inc/myDAT.h"

extern void     ReadHbscFile(char *fileName, data_t * data);
extern void     SetHbscData(int num1, char *atom1, int num2, char *atom2, double y, double sigma, data_t *data);

extern double   CalcChi2_Hbsc_Fwd(polyPeptide_t * polyP, data_t * data, int ppn);
extern double   CalcChi2_Hbsc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, int ppn);
extern double   CalcChi2_Hbsc_Rev(polyPeptide_t * polyP, data_t * data, int ppn);
extern double   CalcChi2_Hbsc_Rev_Flat(polyPeptide_t * polyP, data_t * data, int ppn);

extern double   DisplayChi2_Hbsc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, int ppn);
extern double   DisplayChi2_Hbsc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, int ppn);

#endif
