#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include "../inc/myPARAMS.h"
#include "../inc/mySTRUCT.h"
#include "../inc/myGEOMETRY.h"
#include "../inc/myRDC.h"
#include "../inc/myDAT.h"
#include "../inc/myRAMACHANDRAN.h"
#include "../inc/myPPGEO.h"
#include "../inc/myHBSC.h"
#include "../inc/myMINIMISATION.h"
#include "../inc/gradient.h"
#include "../inc/myJUNC.h"

#include "../inc/meccano.h"

#define PARAM_NB (3 * pPlaneNbInMin)

/******************************
 *    Structure Definitions   *
 ******************************/

typedef struct {
    gsl_vector     *x;
    double          chisq;
} variable_t;


/******************************
 *    Global Variables        *
 ******************************/

int             printFlag = 0;

int             maxValueBest;
int             minValueBest;
int             minValues;
int             ppNbMin;

/*
 * Variables for Random Distributions 
 */
const gsl_rng_type *type;
gsl_rng        *r;

/******************************
 *    Function declaration    *
 ******************************/

void            DelParams(param_t * params);
void            InitRandom(void);
void            DataInit(char **media, int mediaNb, param_t * params);
void            InitExtras(int firstPPlaneFrag, int lastPPlaneFrag, param_t *params);
int             MediaReading(char *mediaFileName, char ***tensorStatFileNames, char ***tensorDynFileNames, char ***rdcFileNames);
void            LectureFichiers(char **tensorDynFileNames, char **tensorStatFileNames, char **media, char *csiFileName, char *talosFileName, char *hbscFileName, char *sequenceFileName, char *ramaDBFileName, param_t * params);
double          Rotation_f(const gsl_vector * x, void *params);
void            Rotation_df(const gsl_vector * x, void *params, gsl_vector * df);
void            Rotation_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df);
void            DisplayChi2(const gsl_vector * x, void *params, FILE * chi2OutputFile);
int             CompareChi2(const void *a, const void *b);
int             TestBest(int nn, variable_t Best[]);
double          DeepMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params);
double          FastMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params);
void            PeptideChainConstruction(param_t *params);
void            WriteOutputFile(char *outputFileName, param_t *params);

void RunMeccanoFromFiles(char *outputFileName, char *mediaFileName,
    char *csiFileName, char *talosFileName, char *hbscFileName,
    char *sequenceFileName, char *ramaDBFileName, int firstPPlaneFrag,
    int lastPPlaneFrag, int ppNbMin_, int minValueBest_, int maxValueBest_)
{
    char          **tensorDynFileNames;
    char          **tensorStatFileNames;
    int             mediaNb;
    char          **media;
    param_t        *params = (param_t *) malloc(sizeof(param_t));

    ppNbMin = ppNbMin_;
    minValueBest = minValueBest_;
    maxValueBest = maxValueBest_;

    minValues = 5 * minValueBest;

    InitRandom();

    mediaNb = MediaReading(mediaFileName, &tensorDynFileNames, &tensorStatFileNames, &media);

    DataInit(media, mediaNb, params);

    LectureFichiers(tensorDynFileNames, tensorStatFileNames, media, csiFileName, talosFileName, hbscFileName, sequenceFileName, ramaDBFileName, params);
    
    if (lastPPlaneFrag - firstPPlaneFrag < 0) {
        fprintf(stderr, "\n Error: lastPPlaneFrag < firstPPlaneFrag\n\n");
        exit(0);
    }
    
    InitExtras(firstPPlaneFrag, lastPPlaneFrag, params);
    
    fprintf(stdout, "\n=====================================================================================\n");
    PeptideChainConstruction(params);
    fprintf(stdout, "\n=====================================================================================\n");

    WriteOutputFile(outputFileName, params);
    WritePdbFile(outputFileName, params->polyP);

    DelParams(params);
}

void RunMeccanoFromParams(int firstPPlaneFrag, int lastPPlaneFrag, int ppNbMin_,
        int minValueBest_, int maxValueBest_, param_t *params)
{
    ppNbMin = ppNbMin_;
    minValueBest = minValueBest_;
    maxValueBest = maxValueBest_;

    minValues = 5 * minValueBest;

    InitRandom();

    InitExtras(firstPPlaneFrag, lastPPlaneFrag, params);

    PeptideChainConstruction(params);
}

param_t *NewParams()
{
    param_t        *params = (param_t *) malloc(sizeof(param_t));

    params->pPlaneRef = NULL;
    params->tensorStat = NULL;
    params->tensorDyn = NULL;
    params->data = NULL;
    params->polyP = NULL;
    params->ramaDB = NULL;

    return params;
}

void DelParams(param_t * params)
{
    if (params->pPlaneRef)
        DelPepPla(params->pPlaneRef);

    if (params->polyP)
        DelPolyPeptide(params->polyP);

    if (params->ramaDB)
        DelRamaDB(params->ramaDB);

    if (params->tensorStat && params->data)
        DelTensor(params->tensorStat, params->data->mediaNb);

    if (params->tensorDyn && params->data)
        DelTensor(params->tensorDyn, params->data->mediaNb);

    if (params->data)
        DelData(params->data);

    free(params);
}

void InitRandom(void)
/*
 * Initialization of "random" variables 
 */
{
    time_t          t1;

    (void) time(&t1);
    type = gsl_rng_default;
    r = gsl_rng_alloc(type);
    gsl_rng_set(r, (unsigned long int) t1);
}

void DataInit(char **media, int mediaNb, param_t * params)

/*
 * Initialization of data_t and tensor_t
 */
{
    int             totalPepPlaNb, firstPepPla;

    TotalPepPlaNbDet(media, mediaNb, &totalPepPlaNb, &firstPepPla);

    params->pPlaneRef = InitPepPlaFwd();

    params->polyP = InitPolyPeptide(totalPepPlaNb, firstPepPla);

    params->data = InitData(mediaNb, totalPepPlaNb, firstPepPla);

    params->tensorStat = InitTensor(mediaNb);
    params->tensorDyn = InitTensor(mediaNb);
}

void InitExtras(int firstPPlaneFrag, int lastPPlaneFrag, param_t *params)
{
    if (firstPPlaneFrag < params->polyP->firstPepPlaNb) {
        fprintf(stderr, "\n Warning: firstPPlaneFrag < firstPPlane\n\n");
        firstPPlaneFrag = params->polyP->firstPepPlaNb;
    }
    if (firstPPlaneFrag > (params->data->totalPepPlaNb + params->polyP->firstPepPlaNb)) {
        fprintf(stderr, "\n Error: firstPPlaneFrag > lastPPlane\n\n");
        exit(0);
    }
    if (lastPPlaneFrag > (params->data->totalPepPlaNb + params->polyP->firstPepPlaNb - 1)) {
        fprintf(stderr, "\n Warning: lastPPlaneFrag > lastPPlane\n\n");
        lastPPlaneFrag = (params->data->totalPepPlaNb + params->polyP->firstPepPlaNb - 1);
    }

    params->polyP->windowSize = lastPPlaneFrag - firstPPlaneFrag + 1;
    params->polyP->windowOffset = firstPPlaneFrag - params->polyP->firstPepPlaNb;
    params->polyP->tetFlag = 1;
    params->polyP->tet0 = DEG2RAD(111.0);
    params->polyP->dTet0 = DEG2RAD(4.0);
    params->polyP->tetEr = DEG2RAD(4.0);
}

/***********************
 *    File Readding   *
 ***********************/

int MediaReading(char *mediaFileName, char ***tensorDynFileNames, char ***tensorStatFileNames, char ***rdcFileNames) {
    FILE           *file;
    char            fileName1[1000], fileName2[1000], fileName3[1000];
    int             m;

    /* 
     * Counting of the different media in "media.txt" file 
     */

    int             mediaNb = 0;

    if ((file = fopen(mediaFileName, "r")) == NULL) {
        printf("Ouverture du fichier %s impossible.\n", mediaFileName);
        exit(1);
    }

    while (!feof(file)) {//     if ((params->polyP->windowSize) > (params->data->totalPepPlaNb)) {
//         printf("\n Warning: Window > NbVal    /    NbVal = %d\n\n", params->data->totalPepPlaNb);
//         params->polyP->windowSize = params->data->totalPepPlaNb;
//     }
        // 
//     offsetMax = params->data->totalPepPlaNb - params->polyP->windowSize;

        fscanf(file, "%s %s %s\n", fileName1, fileName2, fileName3);
        mediaNb++;
    }

    fclose(file);

    /* 
     * Reading of the different angle files in "angles.txt" file 
     */

    *tensorDynFileNames = (char **) malloc(mediaNb * sizeof(char *));
    *tensorStatFileNames = (char **) malloc(mediaNb * sizeof(char *));
    *rdcFileNames = (char **) malloc(mediaNb * sizeof(char *));

    if ((file = fopen(mediaFileName, "r")) == NULL) {
        printf("Ouverture du fichier %s impossible.\n", "media.txt");
        exit(1);
    }

    for (m = 0; m < mediaNb; m++) {
        fscanf(file, "%s %s %s\n", fileName1, fileName2, fileName3);

        (*tensorDynFileNames)[m] = (char *) malloc((strlen(fileName1) + 1) * sizeof(char));
        (*tensorStatFileNames)[m] = (char *) malloc((strlen(fileName2) + 1) * sizeof(char));
        (*rdcFileNames)[m] = (char *) malloc((strlen(fileName3) + 1) * sizeof(char));

        strcpy((*tensorDynFileNames)[m], fileName1);
        strcpy((*tensorStatFileNames)[m], fileName2);
        strcpy((*rdcFileNames)[m], fileName3);
    }
    fclose(file);

    return mediaNb;
}

void LectureFichiers(char **tensorDynFileNames, char **tensorStatFileNames, char **media, char *csiFileName, char *talosFileName, char *hbscFileName, char *sequenceFileName, char *ramaDBFileName, param_t * params) {
    int             m;

    /* 
     * sequenceFileNameUENCE
     */
    ReadSequenceFile(sequenceFileName, params->polyP);

    /* 
     * ramaDB_t 
     */
    params->ramaDB = ReadRamaDBFile(ramaDBFileName);
    if (!params->ramaDB)
	exit(0);

    /* 
     * CS ISO 
     */
    ReadCsiFile(csiFileName, params->data);

    /* 
     * tensor_ts - RDCs 
     */
    fprintf(stdout, "*******************\n");
    fprintf(stdout, "* Tensor readding *\n");
    fprintf(stdout, "*******************\n");
    for (m = 0; m < params->data->mediaNb; m++) {
        ReadTensorFile(tensorDynFileNames[m], &(params->tensorDyn[m]));
        ReadTensorFile(tensorStatFileNames[m], &(params->tensorStat[m]));
		PrintTensor(stdout, &(params->tensorStat[m]));

        ReadRdcFile(media[m], m, params->data);
    }
    fprintf(stdout, "*******************\n\n");
    /* 
     * HBSC 
     */
    ReadHbscFile(hbscFileName, params->data);

    /* 
     * phi0 psi0
     */
    ReadPhiPsi(talosFileName, params->data);

}

/***********************************************
 *      Function Definition  |  f(params,x)    *
 ***********************************************/

double Rotation_f(const gsl_vector * x, void *params) {
    double          f = 0.;
    int             j;

    gsl_vector     *vector = gsl_vector_alloc(3);

    int             ppi = ((param_t *) params)->ppi;
    int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
    int             windowStart = ((param_t *) params)->windowStart;
    pPlane_t       *pPlaneRef = ((param_t *) params)->pPlaneRef;
    data_t         *data = ((param_t *) params)->data;
    polyPeptide_t  *polyP = ((param_t *) params)->polyP;
    tensor_t       *tensorStat = ((param_t *) params)->tensorStat;

    for (j = 0; j < pPlaneNbInMin; j++) {
        int             ppn = ppi + j;
        int             firstPepPlaFlag = (windowStart && j == 0);
        gsl_vector_const_view xj = gsl_vector_const_subvector(x, 3 * j, 3);

        OrientatePepPlaFwd(&xj.vector, firstPepPlaFlag, ppn, pPlaneRef, polyP);

        f += CalcChi2Stat_PPRdc /*_Flat*/ (polyP, data, tensorStat, ppn);

        if (firstPepPlaFlag == 0) {

            f += CalcChi2Stat_CbHaRdc /*_Flat*/ (polyP, data, tensorStat, ppn - 1);

            if (data->hnHnFlag != 0)
                f += CalcChi2Stat_HHRdc_Fwd /*_Flat*/ (polyP, data, tensorStat, ppn);
            if (data->haHnFlag != 0)
                f += CalcChi2Stat_HaHRdc_Fwd /*_Flat*/ (polyP, data, tensorStat, ppn);
            if (data->cHnFlag != 0)
                f += CalcChi2Stat_CHRdc_Fwd /*_Flat*/ (polyP, data, tensorStat, ppn);
            if (data->cbHnFlag != 0)
                f += CalcChi2Stat_CbHRdc_Fwd /*_Flat*/ (polyP, data, tensorStat, ppn);
            if (data->hbscFlag != 0)
                f += CalcChi2_Hbsc_Fwd /*_Flat*/ (polyP, data, ppn);

            f += CalcChi2_JuncGeo(polyP, data, ppn - 1);
        }
    }

    gsl_vector_free(vector);

    return f;
}

// void Rotation_df(const gsl_vector * x, void *params, gsl_vector * df) {
//     int             j;
//     double          h, xj, fhp, f;
//     int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
//     gsl_vector     *xhp = gsl_vector_alloc(PARAM_NB);
// 
//     f = Rotation_f(x, params);
// 
//     for (j = 0; j < PARAM_NB; j++) {
//         xj = GET(x, j);
// 
//         if (xj)
//             h = 1e-10 * xj;
//         else
//             h = 1e-10;
// 
//         gsl_vector_memcpy(xhp, x);
//         SET(xhp, j, xj + h);
//         fhp = Rotation_f(xhp, params);
// 
//         SET(df, j, (fhp - f) / h);
//     }
//     gsl_vector_free(xhp);
// }

void Rotation_df(const gsl_vector * x, void *params, gsl_vector * df) {
    gradient(&Rotation_f, x, params, df);
}

void Rotation_fdf(const gsl_vector * x, void *params, double *f, gsl_vector * df) {
    *f = Rotation_f(x, params);
    Rotation_df(x, params, df);
}

void DisplayChi2(const gsl_vector * x, void *params, FILE * chi2OutputFile) {
    int             j;

    gsl_vector     *vector = gsl_vector_alloc(3);

    int             ppi = ((param_t *) params)->ppi;
    int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
    int             windowStart = ((param_t *) params)->windowStart;
    pPlane_t       *pPlaneRef = ((param_t *) params)->pPlaneRef;
    data_t         *data = ((param_t *) params)->data;
    polyPeptide_t  *polyP = ((param_t *) params)->polyP;
    tensor_t       *tensorStat = ((param_t *) params)->tensorStat;

    for (j = 0; j < pPlaneNbInMin; j++) {
        double          f = 0.0;
        double          f_tot = 0.0;

        int             ppn = ppi + j;
        int             firstPepPlaFlag = (windowStart && j == 0);
        gsl_vector_const_view xj = gsl_vector_const_subvector(x, 3 * j, 3);

        OrientatePepPlaFwd(&xj.vector, firstPepPlaFlag, ppn, pPlaneRef, polyP);

        if (firstPepPlaFlag == 0) {
            f = 0.0;
            DisplayJunctionParam(chi2OutputFile, polyP, ppn - 1);
            f += DisplayChi2_JuncGeo(chi2OutputFile, polyP, data, ppn - 1);
            f += DisplayChi2_CbHaRdc(chi2OutputFile, polyP, data, tensorStat, ppn - 1);
            fprintf(chi2OutputFile, "\n");
            fprintf(chi2OutputFile, "Chi2 Junction: %8.3lf\n", f);
            fprintf(chi2OutputFile, "\n-------------------------------------------------------------------------------------\n");
            f_tot += f;
        }

        f =0.0;
        fprintf(chi2OutputFile, "\nPeptide Plane n°%d\n\n", data->pepPla[ppn].pepPlaNb);
        f += DisplayChi2_PPRdc(chi2OutputFile, polyP, data, tensorStat, ppn);
        fprintf(chi2OutputFile, "\n");
        fprintf(chi2OutputFile, "Chi2 Peptide Plane: %8.3lf\n", f);
        fprintf(chi2OutputFile, "\n-------------------------------------------------------------------------------------\n");
        f_tot += f;

        if ((firstPepPlaFlag == 0) && (data->hnHnFlag + data->haHnFlag + data->cHnFlag + data->cbHnFlag + data->hbscFlag != 0)) {
            f =0.0;
            fprintf(chi2OutputFile, "\nLong Range Couplings:\n\n");
            if (data->hnHnFlag != 0)
                f += DisplayChi2_HHRdc_Fwd(chi2OutputFile, polyP, data, tensorStat, ppn);
            if (data->haHnFlag != 0 /* && ppn-1>0 */ )
                f += DisplayChi2_HaHRdc_Fwd(chi2OutputFile, polyP, data, tensorStat, ppn);
            if (data->cHnFlag != 0 /* && ppn-1>0 */ )
                f += DisplayChi2_CHRdc_Fwd(chi2OutputFile, polyP, data, tensorStat, ppn);
            if (data->cbHnFlag != 0 /* && ppn-1>0 */ )
                f += DisplayChi2_CbHRdc_Fwd(chi2OutputFile, polyP, data, tensorStat, ppn);

            if (data->hbscFlag != 0)
                f += DisplayChi2_Hbsc_Fwd(chi2OutputFile, polyP, data, ppn);

            fprintf(chi2OutputFile, "\n");
            fprintf(chi2OutputFile, "Chi2 Long Range: %8.3lf\n", f);
            fprintf(chi2OutputFile, "\n-------------------------------------------------------------------------------------\n");
            f_tot += f;
        }
        fprintf(chi2OutputFile, "\n");
        fprintf(chi2OutputFile, "Chi2 Total: %8.3lf\n", f_tot);
        fprintf(chi2OutputFile, "\n=====================================================================================\n");
    }

    /**************/

    gsl_vector_free(vector);
}

int CompareChi2(const void *a, const void *b)

/*
 * Comparison function for qsort() 
 */
{
    double          aa = ((variable_t *) a)->chisq;
    double          bb = ((variable_t *) b)->chisq;

    if (aa > bb)
        return 1;
    else if (aa < bb)
        return -1;
    else
        return 0;
}

int TestBest(int nn, variable_t Best[])

/*
 * Stop condition of the minimization step 
 */
{
    if (nn < minValues)
        return 1;
    else
        qsort(Best, maxValueBest, sizeof(variable_t), CompareChi2);

    if (Best[0].chisq < 1e-4)
        return 0;
    else if (fabs(Best[0].chisq - Best[minValueBest - 1].chisq) < 1e-4)
        return 0;
    else
        return 1;
}

/*
 * Minimization Step 
 */
double DeepMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params) {
    int             i, j;
    variable_t      Best[maxValueBest];
    int             nn = 0;

    int             ppi = ((param_t *) params)->ppi;
    int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
    int             windowStart = ((param_t *) params)->windowStart;
    polyPeptide_t  *polyP = ((param_t *) params)->polyP;

    double          bestChi2 = 1e16;

    for (j = 0; j < maxValueBest; j++) {
        Best[j].x = gsl_vector_alloc(x->size);
        Best[j].chisq = 1e16;
    }

    nn = 0;
    do {
        printFlag = 0;
        Best[nn].chisq = PSOMinRotFwd(Func_f, Best[nn].x, params);
//         Best[nn].chisq = SimplexMin(Func_f, Best[nn].x, params);
//         Best[nn].chisq = SteepDescMin(Func_f, Func_df, Func_fdf, Best[nn].x, params);
        Best[nn].chisq = ConjGradMin(Func_f, Func_df, Func_fdf, Best[nn].x, params, 50000);
        nn++;
        fprintf(stderr, "|");
    } while (TestBest(nn, Best) && nn < maxValueBest);

    fprintf(stderr, "\n");

    if (Best[0].chisq < bestChi2) {
        bestChi2 = Best[0].chisq;
        gsl_vector_memcpy(x, Best[0].x);
    }

    for (i = 0; i < pPlaneNbInMin; i++) {
        if (windowStart == 0 || i != 0) {
            polyP->acAm[ppi - 1 + i].phi = GET(x, 0 + 3 * i);
            polyP->acAm[ppi - 1 + i].psi = GET(x, 1 + 3 * i);
            polyP->acAm[ppi - 1 + i].tet = GET(x, 2 + 3 * i);
        } else {
            polyP->acAm[ppi - 1].first[0] = GET(x, 0);
            polyP->acAm[ppi - 1].first[1] = GET(x, 1);
            polyP->acAm[ppi - 1].first[2] = GET(x, 2);
        }
    }

    for (j = 0; j < maxValueBest; j++)
        gsl_vector_free(Best[j].x);

    return bestChi2;
}

double FastMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params) {
    int             i;
    double          chi2 = 1e32;

    int             ppi = ((param_t *) params)->ppi;
    int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
    int             windowStart = ((param_t *) params)->windowStart;
    polyPeptide_t  *polyP = ((param_t *) params)->polyP;

    for (i = 0; i < pPlaneNbInMin; i++) {
        if (windowStart == 0 || i != 0) {
            SET(x, 0 + 3 * i, polyP->acAm[ppi - 1 + i].phi);
            SET(x, 1 + 3 * i, polyP->acAm[ppi - 1 + i].psi);
            SET(x, 2 + 3 * i, polyP->acAm[ppi - 1 + i].tet);
        } else {
            SET(x, 0, polyP->acAm[ppi - 1].first[0]);
            SET(x, 1, polyP->acAm[ppi - 1].first[1]);
            SET(x, 2, polyP->acAm[ppi - 1].first[2]);
        }
    }

    printFlag = 1;
//     chi2 = SimplexMin(Func_f, x, params);
//     chi2 = SteepDescMin(Func_f, Func_df, Func_fdf, x, params);
    chi2 = ConjGradMin(Func_f, Func_df, Func_fdf, x, params, 100);
    printFlag = 0;

    for (i = 0; i < pPlaneNbInMin; i++) {
        if (windowStart == 0 || i != 0) {
            polyP->acAm[ppi - 1 + i].phi = GET(x, 0 + 3 * i);
            polyP->acAm[ppi - 1 + i].psi = GET(x, 1 + 3 * i);
            polyP->acAm[ppi - 1 + i].tet = GET(x, 2 + 3 * i);
        } else {
            polyP->acAm[ppi - 1].first[0] = GET(x, 0);
            polyP->acAm[ppi - 1].first[1] = GET(x, 1);
            polyP->acAm[ppi - 1].first[2] = GET(x, 2);
        }
    }

    return chi2;
}


void PeptideChainConstruction(param_t *params) {
    int             j;
    variable_t     *variables = (variable_t *) malloc(sizeof(variable_t));
    
    for (j = 0; j < params->polyP->windowSize; j++) {

        params->ppi = params->polyP->windowOffset + 1 + j;

        if (j == 0) {
            params->windowStart = 1;
            params->pPlaneNbInMin = (ppNbMin > 1 ? ppNbMin : 2);
        } else {
            params->windowStart = 0;
            params->pPlaneNbInMin = ppNbMin;
        }

        if (params->pPlaneNbInMin > params->polyP->windowSize + params->polyP->windowOffset - params->ppi + 1)
            params->pPlaneNbInMin = params->polyP->windowSize + params->polyP->windowOffset - params->ppi + 1;


        /*****************************************/

        variables->x = gsl_vector_alloc(3 * params->pPlaneNbInMin);
        variables->chisq = DeepMin(&Rotation_f, &Rotation_df, &Rotation_fdf, variables->x, params);
        DisplayChi2(variables->x, params, stdout);

        gsl_vector_free(variables->x);

        /*************************************************/
        if ((j + 1) % 10 == 0) {
            params->ppi = params->polyP->windowOffset + 1 + j - 9;
            params->pPlaneNbInMin = 10;
            if (params->ppi == params->polyP->windowOffset + 1)
                params->windowStart = 1;

            variables->x = gsl_vector_alloc(params->pPlaneNbInMin * 3);
            variables->chisq = FastMin(&Rotation_f, &Rotation_df, &Rotation_fdf, variables->x, params);
            gsl_vector_free(variables->x);
        }

//         if ( (j+1)%5 == 0 ) {    
//             params->ppi = params->polyP->windowOffset + 1;
//             params->pPlaneNbInMin = j + 1;
//             params->windowStart = 1;
        // 
//             variables->x = gsl_vector_alloc(params->pPlaneNbInMin * 3);
//             variables->chisq = FastMin(&Rotation_f, &Rotation_df, &Rotation_fdf, variables->x, params);
//             gsl_vector_free(variables->x);
//         }

        /*************************************************/
    }
    params->ppi = params->polyP->windowOffset + 1;
    params->pPlaneNbInMin = params->polyP->windowSize;
    params->windowStart = 1;

    variables->x = gsl_vector_alloc(params->pPlaneNbInMin * 3);
    variables->chisq = FastMin(&Rotation_f, &Rotation_df, &Rotation_fdf, variables->x, params);

    free(variables);
}

void WriteOutputFile(char *outputFileName, param_t *params){
    char            chi2OutputFileName[1000];
    variable_t     *variables = (variable_t *) malloc(sizeof(variable_t));
    FILE           *chi2OutputFile;
    
    params->ppi = params->polyP->windowOffset + 1;
    params->pPlaneNbInMin = params->polyP->windowSize;
    params->windowStart = 1;
    
    variables->x = gsl_vector_alloc(3 * params->pPlaneNbInMin);
    
    sprintf(chi2OutputFileName, "%s_%d_%d.out", outputFileName, params->polyP->acAm[params->ppi - 1].resNb, params->polyP->acAm[params->ppi - 1 + params->pPlaneNbInMin].resNb);

    if ((chi2OutputFile = fopen(chi2OutputFileName, "w")) == NULL) {
        printf("Ouverture du fichier %s impossible.\n", chi2OutputFileName);
        exit(1);
    }

    DisplayChi2(variables->x, params, chi2OutputFile);

    printf("\nChi2: %lf\n", variables->chisq);
    
    fclose(chi2OutputFile);
    free(variables);
}
