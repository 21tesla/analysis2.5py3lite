#include "../inc/myHBSC.h"

void ReadHbscFile(char *fileName, data_t * data) {
    FILE           *fileIn;

    if (!(fileIn = fopen(fileName, "r"))) {
        printf("Can't Open %s\n", fileName);
        exit(0);
    }

    while (!feof(fileIn)) {
        int             num1, num2;
        double          y, sigma;
        char            atom1[4], atom2[4];

        fscanf(fileIn, "%d %s %d %s %lf %lf\n", &num1, atom1, &num2, atom2, &y, &sigma);
        SetHbscData(num1, atom1, num2, atom2, y, sigma, data);
    }

    fclose(fileIn);
}

void SetHbscData(int num1, char *atom1, int num2, char *atom2, double y, double sigma, data_t *data)
{
    int             i;
    int             offset = data->firstPepPlaNb - 1;

    if (!strcmp(atom1, "N") && !strcmp(atom2, "O")) {
        if ((num1 - offset >= 0) && (num1 - offset <= data->totalPepPlaNb) && (num2 + 1 - offset >= 0) && (num2 + 1 - offset <= data->totalPepPlaNb)) {
            int             inc;

            data->hbscFlag = 1;

            i = num1 - offset;
            inc = data->pepPla[i].hbscNO.couplNb;
            SET(data->pepPla[i].hbscNO.y, inc, y);
            SET(data->pepPla[i].hbscNO.sigma, inc, sigma * 0.5);
            ISET(data->pepPla[i].hbscNO.res, inc, num2 + 1 - offset);
            data->pepPla[i].hbscNO.couplNb++;

            i = num2 + 1 - offset;
            inc = data->pepPla[i].hbscNO.couplNb;
            SET(data->pepPla[i].hbscNO.y, inc, y);
            SET(data->pepPla[i].hbscNO.sigma, inc, sigma * 0.5);
            ISET(data->pepPla[i].hbscNO.res, inc, num1 - offset);
            data->pepPla[i].hbscNO.couplNb++;
        }
    }
}

double CalcChi2_Hbsc_Fwd(polyPeptide_t * polyP, data_t * data, int ppn) {
    double          chi2 = 0.0;

    if (data->hbscFlag != 0) {
        lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
        lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn > ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2 += CHI2(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2 += CHI2(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2 += CHI2(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2 += CHI2(R_ho, 1.700, GET(hbscNO->sigma, k));
                }
            }
        }

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn > ppn2) {
                    R_on = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].n_);
                    R_oh = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].h_);

                    if (R_on < 3.000)
                        chi2 += CHI2(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2 += CHI2(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2 += CHI2(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2 += CHI2(R_oh, 1.700, GET(hbscON->sigma, k));
                }
            }
        }
    }
    return chi2;
}

double CalcChi2_Hbsc_Fwd_Flat(polyPeptide_t * polyP, data_t * data, int ppn) {
    double          chi2 = 0.0;

    if (data->hbscFlag != 0) {
        lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
        lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn > ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2 += CHI2_FLAT(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2 += CHI2_FLAT(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2 += CHI2_FLAT(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2 += CHI2_FLAT(R_ho, 1.700, GET(hbscNO->sigma, k));
                }
            }
        }

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn > ppn2) {
                    R_on = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].n_);
                    R_oh = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].h_);

                    if (R_on < 3.000)
                        chi2 += CHI2_FLAT(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2 += CHI2_FLAT(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2 += CHI2_FLAT(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2 += CHI2_FLAT(R_oh, 1.700, GET(hbscON->sigma, k));
                }
            }
        }
    }
    return chi2;
}

double CalcChi2_Hbsc_Rev(polyPeptide_t * polyP, data_t * data, int ppn) {
    double          chi2 = 0.0;

    if (data->hbscFlag != 0) {
        lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
        lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn < ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2 += CHI2(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2 += CHI2(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2 += CHI2(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2 += CHI2(R_ho, 1.700, GET(hbscNO->sigma, k));
                }
            }
        }

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn < ppn2) {
                    R_on = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].n_);
                    R_oh = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].h_);

                    if (R_on < 3.000)
                        chi2 += CHI2(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2 += CHI2(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2 += CHI2(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2 += CHI2(R_oh, 1.700, GET(hbscON->sigma, k));
                }
            }
        }
    }
    return chi2;
}

double CalcChi2_Hbsc_Rev_Flat(polyPeptide_t * polyP, data_t * data, int ppn) {
    double          chi2 = 0.0;

    if (data->hbscFlag != 0) {
        lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
        lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn < ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2 += CHI2_FLAT(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2 += CHI2_FLAT(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2 += CHI2_FLAT(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2 += CHI2_FLAT(R_ho, 1.700, GET(hbscNO->sigma, k));
                }
            }
        }

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn < ppn2) {
                    R_on = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].n_);
                    R_oh = Dist(polyP->acAm[ppn - 1].o_, polyP->acAm[ppn2].h_);

                    if (R_on < 3.000)
                        chi2 += CHI2_FLAT(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2 += CHI2_FLAT(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2 += CHI2_FLAT(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2 += CHI2_FLAT(R_oh, 1.700, GET(hbscON->sigma, k));
                }
            }
        }
    }
    return chi2;
}

double DisplayChi2_Hbsc_Fwd(FILE * output, polyPeptide_t * polyP, data_t * data, int ppn) {
    int             flag = 0;
    double          chi2 = 0.0;
    lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
    lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

    if (hbscNO->couplNb || hbscON->couplNb)
        flag++;

    if (flag != 0) {
        fprintf(output, "%-15s %-20s", "HBSC:", "");

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho, chi2NO = 0.0;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn > ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2NO += CHI2(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2NO += CHI2(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2NO += CHI2(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2NO += CHI2(R_ho, 1.700, GET(hbscNO->sigma, k));

                    chi2 += chi2NO;

                    fprintf(output, "(%-3d %3d)   N-0: %6.1lf  H-O: %6.1lf   Chi2: %6.1lf", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, R_no, R_ho, chi2NO);
                }
            }
        }
        fprintf(output, "\n");

        fprintf(output, "%-15s %-20s", "", "");

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh, chi2ON = 0.0;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn < ppn2) {
                    R_on = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_oh = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_on < 3.000)
                        chi2ON += CHI2(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2ON += CHI2(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2ON += CHI2(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2ON += CHI2(R_oh, 1.700, GET(hbscON->sigma, k));

                    chi2 += chi2ON;

                    fprintf(output, "(%-3d %3d)   O-N: %6.1lf  O-H: %6.1lf   Chi2: %6.1lf", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, R_on, R_oh, chi2ON);
                }
            }
        }
    }
    fprintf(output, "\n\n");
    return chi2;
}

double DisplayChi2_Hbsc_Rev(FILE * output, polyPeptide_t * polyP, data_t * data, int ppn) {
    int             flag = 0;
    double          chi2 = 0.0;
    lgRgCoupl_t    *hbscNO = &(data->pepPla[ppn].hbscNO);
    lgRgCoupl_t    *hbscON = &(data->pepPla[ppn].hbscON);

    if (hbscNO->couplNb || hbscNO->couplNb)
        flag++;

    if (flag != 0) {
        fprintf(output, "%-15s %-20s", "HBSC:", "");

        if (hbscNO->couplNb) {
            int             k;
            double          R_no, R_ho, chi2NO = 0.0;

            for (k = 0; k < hbscNO->couplNb; k++) {
                int             ppn2 = IGET(hbscNO->res, k);

                if (ppn < ppn2) {
                    R_no = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_ho = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_no < 3.000)
                        chi2NO += CHI2(R_no, 3.000, GET(hbscNO->sigma, k));
                    else if (R_no > 3.300)
                        chi2NO += CHI2(R_no, 3.300, GET(hbscNO->sigma, k));

                    if (R_ho > 2.300)
                        chi2NO += CHI2(R_ho, 2.300, GET(hbscNO->sigma, k));
                    else if (R_ho > 1.700)
                        chi2NO += CHI2(R_ho, 1.700, GET(hbscNO->sigma, k));

                    chi2 += chi2NO;

                    fprintf(output, "(%-3d %3d)   N-0: %6.1lf  H-O: %6.1lf   Chi2: %6.1lf", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, R_no, R_ho, chi2NO);
                }
            }
        }
        fprintf(output, "\n");

        fprintf(output, "%-15s %-20s", "", "");

        if (hbscON->couplNb) {
            int             k;
            double          R_on, R_oh, chi2ON = 0.0;

            for (k = 0; k < hbscON->couplNb; k++) {
                int             ppn2 = IGET(hbscON->res, k);

                if (ppn < ppn2) {
                    R_on = Dist(polyP->acAm[ppn].n_, polyP->acAm[ppn2 - 1].o_);
                    R_oh = Dist(polyP->acAm[ppn].h_, polyP->acAm[ppn2 - 1].o_);

                    if (R_on < 3.000)
                        chi2ON += CHI2(R_on, 3.000, GET(hbscON->sigma, k));
                    else if (R_on > 3.300)
                        chi2ON += CHI2(R_on, 3.300, GET(hbscON->sigma, k));

                    if (R_oh > 2.300)
                        chi2ON += CHI2(R_oh, 2.300, GET(hbscON->sigma, k));
                    else if (R_oh > 1.700)
                        chi2ON += CHI2(R_oh, 1.700, GET(hbscON->sigma, k));

                    chi2 += chi2ON;

                    fprintf(output, "(%-3d %3d)   O-N: %6.1lf  O-H: %6.1lf   Chi2: %6.1lf", data->pepPla[ppn].pepPlaNb, data->pepPla[ppn2].pepPlaNb, R_on, R_oh, chi2ON);
                }
            }
        }
    }
    fprintf(output, "\n\n");
    return chi2;
}
