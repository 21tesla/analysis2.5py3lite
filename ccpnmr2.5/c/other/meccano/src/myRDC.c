#include "../inc/myRDC.h"

/************************************************
 * Calculate RDC without cst                    *
 ************************************************/

const double    D_HA = -1.6783854070353251839e-12 * GC * GH * 0.0001 / (R_HA * R_HA * R_HA);
const double    D_CB = -1.6783854070353251839e-12 * GC * GC * 0.0001 / (R_CB * R_CB * R_CB);
const double    D_CC = -1.6783854070353251839e-12 * GC * GC * 0.0001 / (R_CC * R_CC * R_CC);
const double    D_CH = -1.6783854070353251839e-12 * GC * GH * 0.0001 / (R_CH * R_CH * R_CH);
const double    D_CN = -1.6783854070353251839e-12 * GC * GN * 0.0001 / (R_CN * R_CN * R_CN);
const double    D_NH = -1.6783854070353251839e-12 * GN * GH * 0.0001 / (R_NH * R_NH * R_NH);
const double    D_CA1H = -1.6783854070353251839e-12 * GC * GH * 0.0001 / (R_CA1H * R_CA1H * R_CA1H);
const double    D_CA2H = -1.6783854070353251839e-12 * GC * GH * 0.0001 / (R_CA2H * R_CA2H * R_CA2H);
const double    D_HH = -1.6783854070353251839e-12 * GH * GH * 0.0001;
const double    D_C_H = -1.6783854070353251839e-12 * GC * GH * 0.0001;

tensor_t       *InitTensor(int mediaNb) {
    int             m;
    tensor_t       *tensor = (tensor_t *) malloc(mediaNb * sizeof(tensor_t));

    for (m = 0; m < mediaNb; m++) {
        tensor[m].rotMat = gsl_matrix_alloc(3, 3);
        tensor[m].xCol = gsl_vector_alloc(3);
        tensor[m].zCol = gsl_vector_alloc(3);
    }
    strcpy(tensor[m].name, "");

    return tensor;
}

void ReadTensorFile(char *tensorName, tensor_t * tensor) {
    double          Aa, Ar, alpha, beta, gamma;
    gsl_matrix     *matrix;

    char           *ptr;
    FILE           *fileIn;

    if (!(fileIn = fopen(tensorName, "r"))) {
        printf("Can't open %s\n", tensorName);
        exit(0);
    }

    fscanf(fileIn, "%lf %lf %lf %lf %lf  \n", &Aa, &Ar, &alpha, &beta, &gamma);

    tensor->Aa = Aa;
    tensor->Ar = Ar;
    tensor->alpha = DEG2RAD(alpha);
    tensor->beta = DEG2RAD(beta);
    tensor->gamma = DEG2RAD(gamma);

    ptr = strrchr(tensorName, '/');
    if (ptr)
        strcpy(tensor->name, ptr + 1);
    else
        strcpy(tensor->name, tensorName);

    matrix = MatEuler(DEG2RAD(alpha), DEG2RAD(beta), DEG2RAD(gamma));
    gsl_matrix_transpose_memcpy(tensor->rotMat, matrix);
    gsl_matrix_get_col(tensor->xCol, matrix, 0);
    gsl_matrix_get_col(tensor->zCol, matrix, 2);

    fclose(fileIn);
}

void SetTensor(char *tensorName, double Aa, double Ar, double alpha, double beta, double gamma, tensor_t * tensor) {
    gsl_matrix     *matrix;

    strcpy(tensor->name, tensorName);

    tensor->Aa = Aa;
    tensor->Ar = Ar;
    tensor->alpha = alpha;
    tensor->beta = beta;
    tensor->gamma = gamma;

    matrix = MatEuler(alpha, beta, gamma);
    gsl_matrix_transpose_memcpy(tensor->rotMat, matrix);
    gsl_matrix_get_col(tensor->xCol, matrix, 0);
    gsl_matrix_get_col(tensor->zCol, matrix, 2);
}

void PrintTensor(FILE * output, tensor_t * tensor) {
    fprintf(output, "TENSOR: %s\n", tensor->name);
    fprintf(output, "    %-4s %7.3lf    %-6s %8.3lf\n", "Aa:", tensor->Aa, "Alpha:", RAD2DEG(tensor->alpha));
    fprintf(output, "    %-4s %7.3lf    %-6s %8.3lf\n", "Ar:", tensor->Ar, "Beta:", RAD2DEG(tensor->beta));
    fprintf(output, "    %-4s %7s    %-6s %8.3lf\n", "", "", "Gamma:", RAD2DEG(tensor->gamma));
}

void TotalPepPlaNbDet(char **media, int mediaNb, int *totalPepPlaNb, int *firstPepPlaNb) {
    FILE           *fileIn;
    int             m;

    int             num1, num2;
    double          y;
    char            atom1[10], atom2[10], end[200];

    int             lPepPla = 0, fPepPla = 1000000;

    for (m = 0; m < mediaNb; m++) {
        if (!(fileIn = fopen(media[m], "r"))) {
            printf("\n\n%s\n", media[m]);
            exit(0);
        }

        while (!feof(fileIn)) {
            fscanf(fileIn, "%d %s %d %s %lf %[^\n] \n", &num1, atom1, &num2, atom2, &y, end);
	    UpdatePepPlaNbDet(num1, atom1, num2, atom2, &fPepPla, &lPepPla);
        }
        fclose(fileIn);
    }
    (*totalPepPlaNb) = lPepPla - fPepPla + 1;
    (*firstPepPlaNb) = fPepPla;

    printf("\n\nNumber of the first peptide plane: %3d\n", fPepPla);
    printf("Number of the last peptide plane:  %3d\n", lPepPla);
    printf("Total number of peptide planes:    %3d\n\n", (*totalPepPlaNb));
}

void UpdatePepPlaNbDet(int num1, char *atom1, int num2, char *atom2, int *fPepPla, int *lPepPla) {
    if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H"))
      && !strcmp(atom2, "N") && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H"))
      && !strcmp(atom1, "N") && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if (!strcmp(atom1, "C") && !strcmp(atom2, "N")
      && (num2 - num1) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1 + 1);
        *lPepPla = GSL_MAX(*lPepPla, num1 + 1);
    } else if (!strcmp(atom2, "C") && !strcmp(atom1, "N")
      && (num1 - num2) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if (!strcmp(atom1, "C") && !strcmp(atom2, "CA")
      && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1 + 1);
        *lPepPla = GSL_MAX(*lPepPla, num1 + 1);
    } else if (!strcmp(atom2, "C") && !strcmp(atom1, "CA")
      && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1 + 1);
        *lPepPla = GSL_MAX(*lPepPla, num1 + 1);
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H"))
      && !strcmp(atom2, "CA") && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H"))
      && !strcmp(atom1, "CA") && !(num1 - num2)) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H"))
      && !strcmp(atom2, "CA") && (num1 - num2) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H"))
      && !strcmp(atom1, "CA") && (num2 - num1) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1 + 1);
        *lPepPla = GSL_MAX(*lPepPla, num1 + 1);
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H"))
      && !strcmp(atom2, "C") && (num1 - num2) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1);
        *lPepPla = GSL_MAX(*lPepPla, num1);
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H"))
      && !strcmp(atom1, "C") && (num2 - num1) == 1) {
        *fPepPla = GSL_MIN(*fPepPla, num1 + 1);
        *lPepPla = GSL_MAX(*lPepPla, num1 + 1);
    }
}

double RdcStat(gsl_vector * vec, tensor_t * tensor) {
    double          xc, zc;
    double          Dr = tensor->Ar;
    double          Da = tensor->Aa;

    gsl_blas_ddot(vec, tensor->xCol, &xc);
    gsl_blas_ddot(vec, tensor->zCol, &zc);

    return (3.0 * (Da + (Dr / 2.0)) * gsl_pow_2(zc) - (Da + (Dr * 1.5)) + 3.0 * Dr * gsl_pow_2(xc));
}

double RdcDynS(double S, const gsl_vector * vec, const tensor_t * tensor) {
    double          xc, zc;
    double          Dr = tensor->Ar;
    double          Da = tensor->Aa;

    gsl_blas_ddot(vec, tensor->xCol, &xc);
    gsl_blas_ddot(vec, tensor->zCol, &zc);

    return S * (3.0 * (Da + (Dr / 2.0)) * gsl_pow_2(zc) - (Da + (Dr * 1.5)) + 3.0 * Dr * gsl_pow_2(xc));
}

double RdcDynGAF(double sig, const gsl_vector * axis, const gsl_vector * vec, const tensor_t * tensor) {
    double          Dr = tensor->Ar;
    double          Da = tensor->Aa;

    double          theta_axis;
    double          phi_axis;
    double          theta_vec_rot;
    double          phi_vec_rot;
    double          r;

    gsl_vector     *vec_rot = gsl_vector_alloc(3);

    Cart2Spher(axis, &r, &theta_axis, &phi_axis);

    RotEuler(vec, vec_rot, phi_axis, theta_axis, 0.0);

    Cart2Spher(vec_rot, &r, &theta_vec_rot, &phi_vec_rot);

    gsl_vector_free(vec_rot);

    return Gaf1D(Da, Dr, -theta_axis, -phi_axis, theta_vec_rot, -phi_vec_rot, sig);
}

double Gaf1D(double Aa, double Ar, double ta, double pa, double t, double p, double sig) {
    double          sig2 = sig * sig;

    double          s1 = 2. * (3. * cos(t) * cos(t) - 1.);
    double          s2 = 2. * sin(2. * t) * exp(-.5 * sig2);
    double          s3 = 2. * sin(t) * sin(t) * exp(-2. * sig2);

    double          d1 = cos(2. * pa + 2. * p);
    double          d2 = cos(-2. * pa + 2. * p);
    double          d3 = cos(2. * pa + p);
    double          d4 = cos(-2. * pa + p);

    return (0.25 * Aa * (s1 * (3. * gsl_pow_2(cos(ta)) - 1.) + 3. * s2 * sin(2. * ta) * cos(p) + 3. * s3 * gsl_pow_2(sin(ta)) * cos(2. * p)) + 0.375 * Ar * (s1 * gsl_pow_2(sin(ta)) * cos(2. * pa) - 2. * s2 * sin(ta) * (d3 * gsl_pow_2(cos(.5 * ta)) - d4 * gsl_pow_2(sin(.5 * ta))) + 2. * s3 * (d1 * gsl_pow_4(cos(.5 * ta)) + d2 * gsl_pow_4(sin(.5 * ta)))));

}

double CalcChi2Stat_PPRdc(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m;

    gsl_vector     *vector = gsl_vector_alloc(3);

    acAm_t         *acAmPrev = &(polyP->acAm[ppn - 1]);
    acAm_t         *acAm = &(polyP->acAm[ppn]);

    double          csa11 = data->pepPla[ppn].csa11;
    double          csa22 = data->pepPla[ppn].csa22;
    int             csiFlag = data->pepPla[ppn].csiFlag;

    double          chi2 = 0.0;

    for (m = 0; m < data->mediaNb; m++) {
        gsl_vector     *y = data->pepPla[ppn].medium[m].y;
        gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
        gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

        if (IGET(flag, CC_)) {
            Points2UVector(acAmPrev->c_, acAmPrev->ca, vector);
            chi2 += CHI2(D_CC * RdcStat(vector, &(tensor[m])), GET(y, CC_), GET(sigma, CC_));
        }
        if (IGET(flag, CH_)) {
            Points2UVector(acAmPrev->c_, acAm->h_, vector);
            chi2 += CHI2(D_CH * RdcStat(vector, &(tensor[m])), GET(y, CH_), GET(sigma, CH_));
        }
        if (IGET(flag, C1H)) {
            Points2UVector(acAmPrev->ca, acAm->h_, vector);
            chi2 += CHI2(D_CA1H * RdcStat(vector, &(tensor[m])), GET(y, C1H), GET(sigma, C1H));
        }
        if (IGET(flag, C2H)) {
            Points2UVector(acAm->ca, acAm->h_, vector);
            chi2 += CHI2(D_CA2H * RdcStat(vector, &(tensor[m])), GET(y, C2H), GET(sigma, C2H));
        }
        if (IGET(flag, NH_)) {
            Points2UVector(acAm->n_, acAm->h_, vector);
            chi2 += CHI2(D_NH * RdcStat(vector, &(tensor[m])), GET(y, NH_), GET(sigma, NH_));
        }
        if (IGET(flag, CN_)) {
            Points2UVector(acAmPrev->c_, acAm->n_, vector);
            chi2 += CHI2(D_CN * RdcStat(vector, &(tensor[m])), GET(y, CN_), GET(sigma, CN_));
        }

        if (IGET(flag, CSA) && csiFlag) {
            double          csa_cmpa, csa_cmpb, csa;

            csa_cmpa = 0.0001 * (2. * csa11 + csa22) * RdcStat(acAm->zAxis, &(tensor[m])) / 3.0;
            csa_cmpb = 0.0001 * (2. * csa22 + csa11) * RdcStat(acAm->xAxis, &(tensor[m])) / 3.0;

            csa = csa_cmpa + csa_cmpb;
            chi2 += CHI2(csa, GET(y, CSA), GET(sigma, CSA));
        }

    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_PPRdc_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m;

    gsl_vector     *vector = gsl_vector_alloc(3);

    acAm_t         *acAmPrev = &(polyP->acAm[ppn - 1]);
    acAm_t         *acAm = &(polyP->acAm[ppn]);

    double          csa11 = data->pepPla[ppn].csa11;
    double          csa22 = data->pepPla[ppn].csa22;
    int             csiFlag = data->pepPla[ppn].csiFlag;

    double          chi2 = 0.0;

    for (m = 0; m < data->mediaNb; m++) {
        gsl_vector     *y = data->pepPla[ppn].medium[m].y;
        gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
        gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

        if (IGET(flag, CC_)) {
            Points2UVector(acAmPrev->c_, acAmPrev->ca, vector);
            chi2 += CHI2_FLAT(D_CC * RdcStat(vector, &(tensor[m])), GET(y, CC_), GET(sigma, CC_));
        }
        if (IGET(flag, CH_)) {
            Points2UVector(acAmPrev->c_, acAm->h_, vector);
            chi2 += CHI2_FLAT(D_CH * RdcStat(vector, &(tensor[m])), GET(y, CH_), GET(sigma, CH_));
        }
        if (IGET(flag, C1H)) {
            Points2UVector(acAmPrev->ca, acAm->h_, vector);
            chi2 += CHI2_FLAT(D_CA1H * RdcStat(vector, &(tensor[m])), GET(y, C1H), GET(sigma, C1H));
        }
        if (IGET(flag, C2H)) {
            Points2UVector(acAm->ca, acAm->h_, vector);
            chi2 += CHI2_FLAT(D_CA2H * RdcStat(vector, &(tensor[m])), GET(y, C2H), GET(sigma, C2H));
        }
        if (IGET(flag, NH_)) {
            Points2UVector(acAm->n_, acAm->h_, vector);
            chi2 += CHI2_FLAT(D_NH * RdcStat(vector, &(tensor[m])), GET(y, NH_), GET(sigma, NH_));
        }
        if (IGET(flag, CN_)) {
            Points2UVector(acAmPrev->c_, acAm->n_, vector);
            chi2 += CHI2_FLAT(D_CN * RdcStat(vector, &(tensor[m])), GET(y, CN_), GET(sigma, CN_));
        }

        if (IGET(flag, CSA) && csiFlag) {
            double          csa_cmpa, csa_cmpb, csa;

            csa_cmpa = 0.0001 * (2. * csa11 + csa22) * RdcStat(acAm->zAxis, &(tensor[m])) / 3.0;
            csa_cmpb = 0.0001 * (2. * csa22 + csa11) * RdcStat(acAm->xAxis, &(tensor[m])) / 3.0;

            csa = csa_cmpa + csa_cmpb;
            chi2 += CHI2_FLAT(csa, GET(y, CSA), GET(sigma, CSA));
        }

    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_CbHaRdc(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);
    acAm_t         *acAm = &(polyP->acAm[ppn]);
    double          chi2 = 0.0;

    for (m = 0; m < data->mediaNb; m++) {
        gsl_vector     *y = data->pepPla[ppn].medium[m].y;
        gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
        gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

        if (IGET(flag, CB_)) {
            Points2UVector(acAm->ca, acAm->cb, vector);
            chi2 += CHI2(D_CB * RdcStat(vector, &(tensor[m])), GET(y, CB_), GET(sigma, CB_));
        }
        if (IGET(flag, HA_)) {
            Points2UVector(acAm->ca, acAm->ha, vector);
            chi2 += CHI2(D_HA * RdcStat(vector, &(tensor[m])), GET(y, HA_), GET(sigma, HA_));
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_CbHaRdc_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);
    acAm_t         *acAm = &(polyP->acAm[ppn]);
    double          chi2 = 0.0;

    for (m = 0; m < data->mediaNb; m++) {
        gsl_vector     *y = data->pepPla[ppn].medium[m].y;
        gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
        gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

        if (IGET(flag, CB_)) {
            Points2UVector(acAm->ca, acAm->cb, vector);
            chi2 += CHI2_FLAT(D_CB * RdcStat(vector, &(tensor[m])), GET(y, CB_), GET(sigma, CB_));
        }
        if (IGET(flag, HA_)) {
            Points2UVector(acAm->ca, acAm->ha, vector);
            chi2 += CHI2_FLAT(D_HA * RdcStat(vector, &(tensor[m])), GET(y, HA_), GET(sigma, HA_));
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        if (HnHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

//                     if(Rhh>5.) chi2 += CHI2(Rhh,5.,0.01);
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HHRdc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        if (HnHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2_FLAT(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

                    if (Rhh > 7.)
                        chi2 += CHI2_FLAT(Rhh, 7., 0.01);
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HHRdc_Rev(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        if (HnHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn < ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

//                     if(Rhh>5.) chi2 += CHI2(Rhh,5.,0.01);
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HHRdc_Rev_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        if (HnHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn < ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2_FLAT(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

                    if (Rhh > 7.)
                        chi2 += CHI2_FLAT(Rhh, 7., 0.01);
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HaHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn - 1].medium[m].HaHn);

        if (HaHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn >= ppn2) {
                    Points2Vector(polyP->acAm[ppn - 1].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m])), GET(HaHn->y, k), GET(HaHn->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HaHRdc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn - 1].medium[m].HaHn);

        if (HaHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn >= ppn2) {
                    Points2Vector(polyP->acAm[ppn - 1].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2_FLAT(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m])), GET(HaHn->y, k), GET(HaHn->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HaHRdc_Rev(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn].medium[m].HaHn);

        if (HaHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn <= ppn2) {
                    Points2Vector(polyP->acAm[ppn].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m])), GET(HaHn->y, k), GET(HaHn->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_HaHRdc_Rev_Flat(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn].medium[m].HaHn);

        if (HaHn->couplNb) {
            int             k;
            double          Rhh;

            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn <= ppn2) {
                    Points2Vector(polyP->acAm[ppn].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += CHI2_FLAT(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m])), GET(HaHn->y, k), GET(HaHn->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_CHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *CHn = &(data->pepPla[ppn - 1].medium[m].CHn);
        lgRgCoupl_t    *HnC = &(data->pepPla[ppn].medium[m].HnC);

        if (CHn->couplNb) {
            int             k;
            double          Rch;

            for (k = 0; k < CHn->couplNb; k++) {
                int             ppn2 = IGET(CHn->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn - 1].c_, polyP->acAm[ppn2].h_, vector);
                    Rch = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_C_H / gsl_pow_3(Rch) * RdcStat(vector, &(tensor[m]))), fabs(GET(CHn->y, k)), GET(CHn->sigma, k));
                }
            }
        }

        if (HnC->couplNb) {
            int             k;
            double          Rhc;

            for (k = 0; k < HnC->couplNb; k++) {
                int             ppn2 = IGET(HnC->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].c_, vector);
                    Rhc = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_C_H / gsl_pow_3(Rhc) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnC->y, k)), GET(HnC->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double CalcChi2Stat_CbHRdc_Fwd(polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    double          chi2 = 0.0;
    int             m, k;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *CbHn = &(data->pepPla[ppn - 1].medium[m].CbHn);
        lgRgCoupl_t    *HnCb = &(data->pepPla[ppn].medium[m].HnCb);

        if (CbHn->couplNb) {
            double          Rch;

            for (k = 0; k < CbHn->couplNb; k++) {
                int             ppn2 = IGET(CbHn->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn - 1].cb, polyP->acAm[ppn2].h_, vector);
                    Rch = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_C_H / gsl_pow_3(Rch) * RdcStat(vector, &(tensor[m]))), fabs(GET(CbHn->y, k)), GET(CbHn->sigma, k));
                }
            }
        }

        if (HnCb->couplNb) {
            double          Rhc;

            for (k = 0; k < HnCb->couplNb; k++) {
                int             ppn2 = IGET(HnCb->res, k);

                if (ppn > ppn2) {
                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].cb, vector);
                    Rhc = NormNormalize(vector);

                    chi2 += CHI2(fabs(D_C_H / gsl_pow_3(Rhc) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnCb->y, k)), GET(HnCb->sigma, k));
                }
            }
        }
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_PPRdc(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m;

    gsl_vector     *vector = gsl_vector_alloc(3);

    acAm_t         *acAmPrev = &(polyP->acAm[ppn - 1]);
    acAm_t         *acAm = &(polyP->acAm[ppn]);

    double          csa11 = data->pepPla[ppn].csa11;
    double          csa22 = data->pepPla[ppn].csa22;
    int             csiFlag = data->pepPla[ppn].csiFlag;

    double          chi2 = 0.0, value;

    fprintf(output, "%-15s %-20s %6s %6s %6s %6s %6s %6s %6s\n", "PepPla RDCs:", "Medium", "NH", "CN", "CC", "CA1H", "CA2H", "CH", "CSA");
    for (m = 0; m < data->mediaNb; m++) {
        gsl_vector     *y = data->pepPla[ppn].medium[m].y;
        gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
        gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

        fprintf(output, "%-15s %-20s ", "", tensor[m].name);
        if (IGET(flag, NH_)) {
            Points2UVector(acAm->n_, acAm->h_, vector);
            chi2 += value = CHI2(D_NH * RdcStat(vector, &(tensor[m])), GET(y, NH_), GET(sigma, NH_));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");
        if (IGET(flag, CN_)) {
            Points2UVector(acAmPrev->c_, acAm->n_, vector);
            chi2 += value = CHI2(D_CN * RdcStat(vector, &(tensor[m])), GET(y, CN_), GET(sigma, CN_));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");
        if (IGET(flag, CC_)) {
            Points2UVector(acAmPrev->c_, acAmPrev->ca, vector);
            chi2 += value = CHI2(D_CC * RdcStat(vector, &(tensor[m])), GET(y, CC_), GET(sigma, CC_));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");
        if (IGET(flag, C1H)) {
            Points2UVector(acAmPrev->ca, acAm->h_, vector);
            chi2 += value = CHI2(D_CA1H * RdcStat(vector, &(tensor[m])), GET(y, C1H), GET(sigma, C1H));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");
        if (IGET(flag, C2H)) {
            Points2UVector(acAm->ca, acAm->h_, vector);
            chi2 += value = CHI2(D_CA2H * RdcStat(vector, &(tensor[m])), GET(y, C2H), GET(sigma, C2H));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");
        if (IGET(flag, CH_)) {
            Points2UVector(acAmPrev->c_, acAm->h_, vector);
            chi2 += value = CHI2(D_CH * RdcStat(vector, &(tensor[m])), GET(y, CH_), GET(sigma, CH_));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");

        if (IGET(flag, CSA) && csiFlag) {
            double          csa_cmpa, csa_cmpb, csa;

            csa_cmpa = 0.0001 * (2. * csa11 + csa22) * RdcStat(acAm->zAxis, &(tensor[m])) / 3.0;
            csa_cmpb = 0.0001 * (2. * csa22 + csa11) * RdcStat(acAm->xAxis, &(tensor[m])) / 3.0;

            csa = csa_cmpa + csa_cmpb;
            chi2 += value = CHI2(csa, GET(y, CSA), GET(sigma, CSA));
            fprintf(output, "%6.1lf ", value);
        } else
            fprintf(output, "%6s ", "XXX");

        fprintf(output, "\n");
    }
    fprintf(output, "\n");
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_CbHaRdc(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, caHaFlag = 0;
    gsl_vector     *vector = gsl_vector_alloc(3);
    acAm_t         *acAm = &(polyP->acAm[ppn]);
    double          chi2 = 0.0, value;

    for (m = 0; m < data->mediaNb; m++)
        caHaFlag += IGET(data->pepPla[ppn].medium[m].flag, CB_) + IGET(data->pepPla[ppn].medium[m].flag, HA_);

    if (caHaFlag != 0) {
        fprintf(output, "%-15s %-20s %6s %6s\n", "Junction RDCs:", "Medium", "CB", "HA");

        for (m = 0; m < data->mediaNb; m++) {
            gsl_vector     *y = data->pepPla[ppn].medium[m].y;
            gsl_vector     *sigma = data->pepPla[ppn].medium[m].sigma;
            gsl_vector_int *flag = data->pepPla[ppn].medium[m].flag;

            fprintf(output, "%-15s %-20s ", "", tensor[m].name);
            if (IGET(flag, CB_)) {
                Points2UVector(acAm->ca, acAm->cb, vector);
                chi2 += value = CHI2(D_CB * RdcStat(vector, &(tensor[m])), GET(y, CB_), GET(sigma, CB_));
                fprintf(output, "%6.1lf ", value);
            } else
                fprintf(output, "%6s ", "XXX");
            if (IGET(flag, HA_)) {
                Points2UVector(acAm->ca, acAm->ha, vector);
                chi2 += value = CHI2(D_HA * RdcStat(vector, &(tensor[m])), GET(y, HA_), GET(sigma, HA_));
                fprintf(output, "%6.1lf ", value);
            } else
                fprintf(output, "%6s ", "XXX");

            fprintf(output, "\n");
        }
        fprintf(output, "\n");

    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_HHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    gsl_vector     *vector = gsl_vector_alloc(3);
    double          chi2 = 0.0, value;

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        for (k = 0; k < HnHn->couplNb; k++) {
            int             ppn2 = IGET(HnHn->res, k);

            if (ppn > ppn2)
                flag++;
        }
    }

    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "HnHn RDC", "Medium");
        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn > ppn2) {
                    double          Rhh;

                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += value = CHI2(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

                    fprintf(output, " (%-3d %3d)   %6.1lf   ", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value);
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_HHRdc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    gsl_vector     *vector = gsl_vector_alloc(3);
    double          chi2 = 0.0, value;

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

        for (k = 0; k < HnHn->couplNb; k++) {
            int             ppn2 = IGET(HnHn->res, k);

            if (ppn < ppn2)
                flag++;
        }
    }

    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "HnHn RDC", "Medium");
        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *HnHn = &(data->pepPla[ppn].medium[m].HnHn);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            for (k = 0; k < HnHn->couplNb; k++) {
                int             ppn2 = IGET(HnHn->res, k);

                if (ppn < ppn2) {
                    double          Rhh;

                    Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += value = CHI2(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnHn->y, k)), GET(HnHn->sigma, k));

                    fprintf(output, " (%-3d %3d)   %6.1lf   ", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value);
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_HaHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    gsl_vector     *vector = gsl_vector_alloc(3);
    double          chi2 = 0.0, value;

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn - 1].medium[m].HaHn);

        for (k = 0; k < HaHn->couplNb; k++) {
            int             ppn2 = IGET(HaHn->res, k);

            if (ppn >= ppn2)
                flag++;
        }
    }
    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "HaHn RDC:", "Medium");
        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *HaHn = &(data->pepPla[ppn - 1].medium[m].HaHn);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn >= ppn2) {
                    double          Rhh;

                    Points2Vector(polyP->acAm[ppn - 1].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += value = CHI2(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m])), GET(HaHn->y, k), GET(HaHn->sigma, k));

                    fprintf(output, " (%-3d %3d)   %6.1lf   ", data->pepPla[ppn - 1].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value);
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_HaHRdc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    gsl_vector     *vector = gsl_vector_alloc(3);
    double          chi2 = 0.0, value;

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *HaHn = &(data->pepPla[ppn].medium[m].HaHn);

        for (k = 0; k < HaHn->couplNb; k++) {
            int             ppn2 = IGET(HaHn->res, k);

            if (ppn <= ppn2)
                flag++;
        }
    }

    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "HaHn RDC:", "Medium");
        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *HaHn = &(data->pepPla[ppn].medium[m].HaHn);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            for (k = 0; k < HaHn->couplNb; k++) {
                int             ppn2 = IGET(HaHn->res, k);

                if (ppn <= ppn2) {
                    double          Rhh;

                    Points2Vector(polyP->acAm[ppn].ha, polyP->acAm[ppn2].h_, vector);
                    Rhh = NormNormalize(vector);

                    chi2 += value = CHI2(fabs(D_HH / gsl_pow_3(Rhh) * RdcStat(vector, &(tensor[m]))), fabs(GET(HaHn->y, k)), GET(HaHn->sigma, k));

                    fprintf(output, " (%-3d %3d)   %6.1lf   ", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value);
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_CHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    double          chi2 = 0.0, value;
    gsl_vector     *vector = gsl_vector_alloc(3);

    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *CHn = &(data->pepPla[ppn - 1].medium[m].CHn);
        lgRgCoupl_t    *HnC = &(data->pepPla[ppn].medium[m].HnC);

        for (k = 0; k < CHn->couplNb; k++) {
            int             ppn2 = IGET(CHn->res, k);

            if (ppn > ppn2)
                flag++;
        }
        for (k = 0; k < HnC->couplNb; k++) {
            int             ppn2 = IGET(HnC->res, k);

            if (ppn > ppn2)
                flag++;
        }
    }

    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "CHn RDC:", "Medium");

        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *CHn = &(data->pepPla[ppn - 1].medium[m].CHn);
            lgRgCoupl_t    *HnC = &(data->pepPla[ppn].medium[m].HnC);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            if (CHn->couplNb) {
                double          Rch;

                for (k = 0; k < CHn->couplNb; k++) {
                    int             ppn2 = IGET(CHn->res, k);

                    if (ppn > ppn2) {
                        Points2Vector(polyP->acAm[ppn - 1].c_, polyP->acAm[ppn2].h_, vector);
                        Rch = NormNormalize(vector);

                        chi2 += value = CHI2(fabs(D_C_H / gsl_pow_3(Rch) * RdcStat(vector, &(tensor[m]))), fabs(GET(CHn->y, k)), GET(CHn->sigma, k));

                        fprintf(output, " (%-3d %3d)   %6.1lf %6.1lfA   ", data->pepPla[ppn - 1].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value, Rch);
                    }
                }
            }

            if (HnC->couplNb) {
                double          Rhc;

                for (k = 0; k < HnC->couplNb; k++) {
                    int             ppn2 = IGET(HnC->res, k);

                    if (ppn > ppn2) {
                        Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].c_, vector);
                        Rhc = NormNormalize(vector);

                        chi2 += value = CHI2(fabs(D_C_H / gsl_pow_3(Rhc) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnC->y, k)), GET(HnC->sigma, k));

                        fprintf(output, " (%-3d %3d)   %6.1lf %6.1lfA   ", data->pepPla[ppn2 - 1].pepPlaNb, data->pepPla[ppn].pepPlaNb, value, Rhc);
                    }
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

double DisplayChi2_CbHRdc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, tensor_t * tensor, int ppn) {
    int             m, k, flag = 0;
    double          chi2 = 0.0, value;
    gsl_vector     *vector = gsl_vector_alloc(3);

    printf("ppn: %d\n", data->pepPla[ppn].pepPlaNb);
    for (m = 0; m < data->mediaNb; m++) {
        lgRgCoupl_t    *CbHn = &(data->pepPla[ppn - 1].medium[m].CbHn);
        lgRgCoupl_t    *HnCb = &(data->pepPla[ppn].medium[m].HnCb);

        for (k = 0; k < CbHn->couplNb; k++) {
            int             ppn2 = IGET(CbHn->res, k);

            if (ppn > ppn2)
                flag++;
        }
        for (k = 0; k < HnCb->couplNb; k++) {
            int             ppn2 = IGET(HnCb->res, k);

            if (ppn > ppn2)
                flag++;
        }
    }

    if (flag != 0) {
        fprintf(output, "%-15s %-20s\n", "CbHn RDC:", "Medium");

        for (m = 0; m < data->mediaNb; m++) {
            lgRgCoupl_t    *CbHn = &(data->pepPla[ppn - 1].medium[m].CbHn);
            lgRgCoupl_t    *HnCb = &(data->pepPla[ppn].medium[m].HnCb);

            fprintf(output, "%-15s %-20s", "", tensor[m].name);
            if (CbHn->couplNb) {
                double          Rch;

                for (k = 0; k < CbHn->couplNb; k++) {
                    int             ppn2 = IGET(CbHn->res, k);

                    if (ppn > ppn2) {
                        Points2Vector(polyP->acAm[ppn - 1].cb, polyP->acAm[ppn2].h_, vector);
                        Rch = NormNormalize(vector);

                        chi2 += value = CHI2(fabs(D_C_H / gsl_pow_3(Rch) * RdcStat(vector, &(tensor[m]))), fabs(GET(CbHn->y, k)), GET(CbHn->sigma, k));

                        fprintf(output, " (%-3d %3d)   %6.1lf %6.1lfA  * ", data->pepPla[ppn - 1].pepPlaNb, data->pepPla[ppn2].pepPlaNb, value, Rch);
                    }
                }
            }

            if (HnCb->couplNb) {
                double          Rhc;

                for (k = 0; k < HnCb->couplNb; k++) {
                    int             ppn2 = IGET(HnCb->res, k);

                    if (ppn > ppn2) {
                        Points2Vector(polyP->acAm[ppn].h_, polyP->acAm[ppn2].cb, vector);
                        Rhc = NormNormalize(vector);

                        chi2 += value = CHI2(fabs(D_C_H / gsl_pow_3(Rhc) * RdcStat(vector, &(tensor[m]))), fabs(GET(HnCb->y, k)), GET(HnCb->sigma, k));

                        fprintf(output, " (%-3d %3d)   %6.1lf %6.1lfA  # ", data->pepPla[ppn2].pepPlaNb, data->pepPla[ppn].pepPlaNb, value, Rhc);
                    }
                }
            }
            fprintf(output, "\n");
        }
        fprintf(output, "\n");
    }
    gsl_vector_free(vector);
    return chi2;
}

void DelTensor(tensor_t *tensor, int mediaNb) {
    int             m;

    for (m = 0; m < mediaNb; m++) {
        gsl_matrix_free(tensor[m].rotMat);
        gsl_vector_free(tensor[m].xCol);
        gsl_vector_free(tensor[m].zCol);
    }

    free(tensor);
}

