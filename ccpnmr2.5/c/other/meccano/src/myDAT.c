#include "../inc/myDAT.h"

data_t         *InitData(int mediaNb, int totalPepPlaNb, int firstPepPlaNb) {
    int             i, m;
    data_t         *data = (data_t *) malloc(sizeof(data_t));

    data->firstPepPlaNb = firstPepPlaNb;
    data->totalPepPlaNb = totalPepPlaNb;

    data->mediaNb = mediaNb;

    data->hnHnFlag = 0;
    data->haHnFlag = 0;
    data->cHnFlag = 0;
    data->hnCFlag = 0;
    data->cbHnFlag = 0;
    data->hnCbFlag = 0;
    data->hbscFlag = 0;
    data->phiPsiFlag = 0;

    data->pepPla = (pepPlaData_t *) malloc((totalPepPlaNb + 1) * sizeof(pepPlaData_t));

    for (i = 0; i <= totalPepPlaNb; i++) {
        data->pepPla[i].csiFlag = 0;
        data->pepPla[i].totalCouplNb = 0;
        data->pepPla[i].medium = (medium_t *) malloc((mediaNb) * sizeof(medium_t));

        for (m = 0; m < mediaNb; m++) {
            data->pepPla[i].medium[m].y = gsl_vector_alloc(PP_COUPL_NB);
            data->pepPla[i].medium[m].sigma = gsl_vector_alloc(PP_COUPL_NB);
            data->pepPla[i].medium[m].flag = gsl_vector_int_alloc(PP_COUPL_NB);
            gsl_vector_set_zero(data->pepPla[i].medium[m].y);
            gsl_vector_set_all(data->pepPla[i].medium[m].sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].flag);

            data->pepPla[i].medium[m].HnHn.couplNb = 0;
            data->pepPla[i].medium[m].HnHn.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnHn.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnHn.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].HnHn.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].HnHn.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].HnHn.res);

            data->pepPla[i].medium[m].HaHn.couplNb = 0;
            data->pepPla[i].medium[m].HaHn.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HaHn.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HaHn.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].HaHn.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].HaHn.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].HaHn.res);

            data->pepPla[i].medium[m].CHn.couplNb = 0;
            data->pepPla[i].medium[m].CHn.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].CHn.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].CHn.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].CHn.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].CHn.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].CHn.res);

            data->pepPla[i].medium[m].HnC.couplNb = 0;
            data->pepPla[i].medium[m].HnC.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnC.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnC.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].HnC.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].HnC.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].HnC.res);

            data->pepPla[i].medium[m].CbHn.couplNb = 0;
            data->pepPla[i].medium[m].CbHn.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].CbHn.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].CbHn.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].CbHn.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].CbHn.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].CbHn.res);

            data->pepPla[i].medium[m].HnCb.couplNb = 0;
            data->pepPla[i].medium[m].HnCb.y = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnCb.sigma = gsl_vector_alloc(totalPepPlaNb);
            data->pepPla[i].medium[m].HnCb.res = gsl_vector_int_alloc(totalPepPlaNb);
            gsl_vector_set_zero(data->pepPla[i].medium[m].HnCb.y);
            gsl_vector_set_all(data->pepPla[i].medium[m].HnCb.sigma, 1e32);
            gsl_vector_int_set_zero(data->pepPla[i].medium[m].HnCb.res);
        }

        data->pepPla[i].hbscNO.couplNb = 0;
        data->pepPla[i].hbscNO.y = gsl_vector_alloc(totalPepPlaNb);
        data->pepPla[i].hbscNO.sigma = gsl_vector_alloc(totalPepPlaNb);
        data->pepPla[i].hbscNO.res = gsl_vector_int_alloc(totalPepPlaNb);
        gsl_vector_set_zero(data->pepPla[i].hbscNO.y);
        gsl_vector_set_all(data->pepPla[i].hbscNO.sigma, 1e32);
        gsl_vector_int_set_zero(data->pepPla[i].hbscNO.res);

        data->pepPla[i].hbscNO.couplNb = 0;
        data->pepPla[i].hbscNO.y = gsl_vector_alloc(totalPepPlaNb);
        data->pepPla[i].hbscNO.sigma = gsl_vector_alloc(totalPepPlaNb);
        data->pepPla[i].hbscNO.res = gsl_vector_int_alloc(totalPepPlaNb);
        gsl_vector_set_zero(data->pepPla[i].hbscNO.y);
        gsl_vector_set_all(data->pepPla[i].hbscNO.sigma, 1e32);
        gsl_vector_int_set_zero(data->pepPla[i].hbscNO.res);

        data->pepPla[i].phiFlag = 0;
        data->pepPla[i].phi0 = data->pepPla[i].psi0 = 0.;
        data->pepPla[i].phiEr = data->pepPla[i].psiEr = 1e32;

        data->pepPla[i].psiFlag = 0;
        data->pepPla[i].psi0 = data->pepPla[i].psi0 = 0.;
        data->pepPla[i].psiEr = data->pepPla[i].psiEr = 1e32;

        data->pepPla[i].pepPlaNb = i + firstPepPlaNb - 1;
    }

    return data;
}

void ReadPhiPsi(char *phiPsiFileName, data_t * data) {
    int             num;
    double          r2, r3, r4, r5;
    FILE           *fileIn;

    if (!(fileIn = fopen(phiPsiFileName, "r"))) {
        printf("Can't Open %s", phiPsiFileName);
        exit(0);
    }

    num = -1;
    while (!feof(fileIn)) {
        fscanf(fileIn, "%d %lf %lf %lf %lf \n", &num, &r2, &r3, &r4, &r5);
	SetPhiPsi(num, r2, r3, r4, r5, data);
    }
    fclose(fileIn);
}

void SetPhiPsi(int num, double r2, double r3, double r4, double r5, data_t *data) {
    int             i, offset = data->firstPepPlaNb - 1;

    if ((num - offset < data->totalPepPlaNb) && (num - offset >= 0)
      && (num >= 0)) {
        i = num - offset;
        data->pepPla[i].phi0 = r2;
        data->pepPla[i].phiEr = r4;
        data->pepPla[i].phiFlag = 1;
        data->pepPla[i].psi0 = r3;
        data->pepPla[i].psiEr = r5;
        data->pepPla[i].psiFlag = 1;
    }
}

void ReadRdcFile(char *rdcFileName, int m, data_t * data) {
    int             num1, num2;
    double          y, sigma, S;
    char            atom1[4], atom2[4];

    FILE           *fileIn;

    if (!(fileIn = fopen(rdcFileName, "r"))) {
        printf("%s", rdcFileName);
        exit(0);
    }

    while (!feof(fileIn)) {
        fscanf(fileIn, "%d %s %d %s %lf %lf %lf\n", &num1, atom1, &num2, atom2, &y, &sigma, &S);
        SetRdcData(m, num1, atom1, num2, atom2, y, sigma, data);
    }

    printf("\n");

    fclose(fileIn);
}

void SetRdcData(int m, int num1, char *atom1, int num2, char *atom2, double y, double sigma, data_t *data) {
    int             i = 0, offset = data->firstPepPlaNb - 1;

    if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && !strcmp(atom2, "N") && !(num1 - num2)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, NH_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, NH_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, NH_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && !strcmp(atom1, "N") && !(num1 - num2)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, NH_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, NH_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, NH_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if (!strcmp(atom1, "C") && !strcmp(atom2, "N") && (num2 - num1) == 1) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CN_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CN_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CN_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if (!strcmp(atom2, "C") && !strcmp(atom1, "N") && (num1 - num2) == 1) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CN_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CN_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CN_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if (!strcmp(atom1, "C") && !strcmp(atom2, "CA") && !(num1 - num2)) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CC_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CC_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CC_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if (!strcmp(atom2, "C") && !strcmp(atom1, "CA") && !(num1 - num2)) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CC_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CC_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CC_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && !strcmp(atom2, "CA") && !(num1 - num2)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, C2H, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, C2H, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, C2H, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && !strcmp(atom1, "CA") && !(num1 - num2)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, C2H, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, C2H, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, C2H, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && !strcmp(atom2, "CA") && (num1 - num2) == 1) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, C1H, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, C1H, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, C1H, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && !strcmp(atom1, "CA") && (num2 - num1) == 1) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, C1H, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, C1H, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, C1H, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && !strcmp(atom2, "C") && (num1 - num2) == 1) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CH_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CH_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CH_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && !strcmp(atom1, "C") && (num2 - num1) == 1) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CH_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CH_, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CH_, 1);
            data->pepPla[i].totalCouplNb++;
        }
    } else if ((!strcmp(atom1, "CA") && !strcmp(atom2, "CB") && !(num1 - num2)) && (num1 - offset < data->totalPepPlaNb) && (num1 - offset >= 0)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CB_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CB_, sigma * sqrt(2.0));
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CB_, 1);
        }
    } else if ((!strcmp(atom2, "CA") && !strcmp(atom1, "CB") && !(num1 - num2)) && (num1 - offset < data->totalPepPlaNb) && (num1 - offset >= 0)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CB_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CB_, sigma * sqrt(2.0));
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CB_, 1);
        }
    } else if ((!strcmp(atom1, "CA") && !strcmp(atom2, "HA") && !(num1 - num2)) && (num1 - offset < data->totalPepPlaNb) && (num1 - offset >= 0)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, HA_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, HA_, sigma * sqrt(2.0));
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, HA_, 1);
        }

    } else if ((!strcmp(atom2, "CA") && !strcmp(atom1, "HA") && !(num1 - num2)) && (num1 - offset < data->totalPepPlaNb) && (num1 - offset >= 0)) {
        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, HA_, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, HA_, sigma * sqrt(2.0));
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, HA_, 1);
        }
    } else if ((!strcmp(atom1, "CSA"))) {
        i = num1 + 1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            gsl_vector_set(data->pepPla[i].medium[m].y, CSA, y);
            gsl_vector_set(data->pepPla[i].medium[m].sigma, CSA, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].flag, CSA, 1);
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && (!strcmp(atom2, "HN") || !strcmp(atom2, "H"))) {
        int             inc;

        i = GSL_MIN(num1, num2) - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnHn.res, inc, GSL_MAX(num1, num2) - offset);
            data->pepPla[i].medium[m].HnHn.couplNb++;
            data->hnHnFlag = 1;
        }

        i = GSL_MAX(num1, num2) - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnHn.res, inc, GSL_MIN(num1, num2) - offset);
            data->pepPla[i].medium[m].HnHn.couplNb++;
            data->hnHnFlag = 1;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && (!strcmp(atom2, "HA"))) {
        int             inc;

        i = num2 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HaHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HaHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HaHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HaHn.res, inc, num1 - offset);
            data->pepPla[i].medium[m].HaHn.couplNb++;
            data->haHnFlag = 1;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && (!strcmp(atom1, "HA"))) {
        int             inc;

        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HaHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HaHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HaHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HaHn.res, inc, num2 - offset);
            data->pepPla[i].medium[m].HaHn.couplNb++;
            data->haHnFlag = 1;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && (!strcmp(atom2, "C")) && (num1 - num2) != 1) {
        int             inc;

        i = num2 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].CHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].CHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].CHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].CHn.res, inc, num1 - offset);
            data->pepPla[i].medium[m].CHn.couplNb++;
            data->cHnFlag = 1;
        }

        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnC.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnC.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnC.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnC.res, inc, num2 - offset);
            data->pepPla[i].medium[m].HnC.couplNb++;
            data->hnCFlag = 1;
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && (!strcmp(atom1, "C")) && (num2 - num1) != 1) {
        int             inc;

        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].CHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].CHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].CHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].CHn.res, inc, num2 - offset);
            data->pepPla[i].medium[m].CHn.couplNb++;
            data->cHnFlag = 1;
        }

        i = num2 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnC.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnC.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnC.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnC.res, inc, num2 - offset);
            data->pepPla[i].medium[m].HnC.couplNb++;
            data->hnCFlag = 1;
        }
    } else if ((!strcmp(atom1, "HN") || !strcmp(atom1, "H")) && (!strcmp(atom2, "CB"))) {
    int             inc;

        i = num2 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].CbHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].CbHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].CbHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].CbHn.res, inc, num1 - offset);
            data->pepPla[i].medium[m].CbHn.couplNb++;
            data->cbHnFlag = 1;
            printf(" %4d H %4d CB %8.3lf\n", data->pepPla[num2 - offset].pepPlaNb, data->pepPla[num1 - offset].pepPlaNb, y);
        }

        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnCb.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnCb.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnCb.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnCb.res, inc, num2 - offset);
            data->pepPla[i].medium[m].HnCb.couplNb++;
            data->hnCbFlag = 1;
            printf(" %4d H %4d CB %8.3lf\n", data->pepPla[num1 - offset].pepPlaNb, data->pepPla[num2 - offset].pepPlaNb, y);
        }
    } else if ((!strcmp(atom2, "HN") || !strcmp(atom2, "H")) && (!strcmp(atom1, "CB"))) {
        int             inc;

        i = num1 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].CbHn.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].CbHn.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].CbHn.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].CbHn.res, inc, num2 - offset);
            data->pepPla[i].medium[m].CbHn.couplNb++;
            data->cbHnFlag = 1;
        }

        i = num2 - offset;
        if (i >= 0 && i <= data->totalPepPlaNb) {
            inc = data->pepPla[i].medium[m].HnCb.couplNb;
            gsl_vector_set(data->pepPla[i].medium[m].HnCb.y, inc, y);
            gsl_vector_set(data->pepPla[i].medium[m].HnCb.sigma, inc, sigma);
            gsl_vector_int_set(data->pepPla[i].medium[m].HnCb.res, inc, num1 - offset);
            data->pepPla[i].medium[m].HnCb.couplNb++;
            data->hnCbFlag = 1;
        }
    }
}

void ReadCsiFile(char *csiFileName, data_t * data) {
    int             num;
    double          r;
    int             offset = data->firstPepPlaNb - 1;
    FILE           *fileIn;

    if (!(fileIn = fopen(csiFileName, "r"))) {
        printf("Can't open %s", csiFileName);
        exit(0);
    }

    while (!feof(fileIn)) {
        int             i;

        fscanf(fileIn, "%d %lf\n", &num, &r);

        i = num + 1 - offset;
        if (i <= data->totalPepPlaNb && i >= 0) {
            data->pepPla[i].csiFlag = 1;
            data->pepPla[i].csa11 = 1000.0 * (247.0 - r);
            data->pepPla[i].csa22 = 1000.0 * (2.0 * r - 332.0);
        }
    }
    fclose(fileIn);
}

void DelData(data_t * data) {
    int i, m;
    int totalPepPlaNb = data->totalPepPlaNb;
    int mediaNb = data->mediaNb;

    for (i = 0; i <= totalPepPlaNb; i++) {
        gsl_vector_free(data->pepPla[i].hbscNO.y);
        gsl_vector_free(data->pepPla[i].hbscNO.sigma);
        gsl_vector_int_free(data->pepPla[i].hbscNO.res);

        for (m = 0; m < mediaNb; m++) {
            gsl_vector_free(data->pepPla[i].medium[m].y);
            gsl_vector_free(data->pepPla[i].medium[m].sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].flag);

            gsl_vector_free(data->pepPla[i].medium[m].HnHn.y);
            gsl_vector_free(data->pepPla[i].medium[m].HnHn.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].HnHn.res);

            gsl_vector_free(data->pepPla[i].medium[m].HaHn.y);
            gsl_vector_free(data->pepPla[i].medium[m].HaHn.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].HaHn.res);

            gsl_vector_free(data->pepPla[i].medium[m].CHn.y);
            gsl_vector_free(data->pepPla[i].medium[m].CHn.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].CHn.res);

            gsl_vector_free(data->pepPla[i].medium[m].HnC.y);
            gsl_vector_free(data->pepPla[i].medium[m].HnC.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].HnC.res);

            gsl_vector_free(data->pepPla[i].medium[m].CbHn.y);
            gsl_vector_free(data->pepPla[i].medium[m].CbHn.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].CbHn.res);

            gsl_vector_free(data->pepPla[i].medium[m].HnCb.y);
            gsl_vector_free(data->pepPla[i].medium[m].HnCb.sigma);
            gsl_vector_int_free(data->pepPla[i].medium[m].HnCb.res);
        }

        free(data->pepPla[i].medium);
    }

    free(data->pepPla);
    free(data);
}
