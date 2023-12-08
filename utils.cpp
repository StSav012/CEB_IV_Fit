#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>

#include "cfoo.h"
#include "coff.h"
#include "mins.h"
#include "utils.h"

long GetExp(char *fname, int iNumCols, int iColNum, double *&Iexp, double *&Vexp, bool bRemOffset)
{ /* Gets experimental values 'Iexp' and 'Vexp' from file 'fname',
	checking 'bRemOffset' if there is need to remove offsets*/

    FILE *in = fopen(fname, "r");
    if (!in) {
        printf("Error! No such filename exists!");
        getchar();
        exit(0);
    } else
        printf("Fitting for %s started! \n", fname);
    FILE *params = fopen("fitparameters_new.txt", "a");
    fprintf(params, "Filename: %s \n", fname);
    fclose(params);
    char c = 0;
    long count = 0;

    double tmp;

    while (fscanf(in, "%lf %lf", &tmp, &tmp) != EOF)
        count++;
    freopen(fname, "r", in);
    count = count / iColNum;

    if (Iexp == NULL)
        Iexp = new double[count];
    if (Vexp == NULL)
        Vexp = new double[count];

    for (int i = 0; i < count; i++) {
        for (int j = 0; j < iNumCols; j++)
            if (j == iColNum - 1) {
                fscanf(in, "%lf", &Vexp[i]);
                fscanf(in, "%lf", &Iexp[i]);
            } else {
                fscanf(in, "%lf", &tmp);
                fscanf(in, "%lf", &tmp);
            }
    }
    if (bRemOffset) {
        double dOffset = 0;
        double dLBound = -0.0005;
        double dRBound = 0.0005;
        COff Off(Iexp, Vexp, count);
        Golden broff;
        broff.ax = dLBound;
        broff.bx = dOffset;
        broff.cx = dRBound;
        dOffset = broff.minimize(Off);
        for (int i = 0; i < count; i++)
            Vexp[i] -= dOffset;
        FILE *ivoff = fopen("IVoffset.txt", "w+");
        for (int i = 0; i < count; i++)
            fprintf(ivoff, "%e %e\n", Vexp[i], Iexp[i]);
        fclose(ivoff);
    }

    return count;
}

void Resample(long countnum,
              long countexp,
              double *&Irex,
              double *&Vrex,
              double *&Iexp,
              double *&Vexp,
              double *Inum,
              double *Vnum)
{ /*Adaptation of experimental data for the fitting algorithm, resamples IV-curve for the specified range*/
    long lLastPos = 0;
    long lCurPos = 0, lLeft = 0, lRight = 0;
    double dCurVal = 0;
    double tmp = 0, tmp1 = 0;
    for (long i = 0; i < countnum; i++) {
        dCurVal = Vnum[i];
        Vrex[i] = Vnum[i];
        lCurPos = 0;
        tmp = std::abs(dCurVal - Vexp[0]);
        for (long j = 1; j < countexp; j++) {
            tmp1 = std::abs(dCurVal - Vexp[j]);
            if (tmp1 < tmp) {
                tmp = tmp1;
                lCurPos = j;
            }
        }
        if ((dCurVal - Vexp[lCurPos]) * (dCurVal - Vexp[lCurPos + 1]) <= 0
            && lCurPos < countexp - 1) {
            lLeft = lCurPos;
            lRight = lCurPos + 1;
        } else {
            lLeft = lCurPos - 1;
            lRight = lCurPos;
        }
        if (lLeft < 0)
            Irex[i] = Iexp[0];
        else if (lRight > countexp - 1)
            Irex[i] = Iexp[countexp - 1];
        else
            Irex[i] = Iexp[lLeft]
                      + (dCurVal - Vexp[lLeft]) * (Iexp[lRight] - Iexp[lLeft])
                            / (Vexp[lRight] - Vexp[lLeft]);
    }
    FILE *outres = fopen("IVresampled.txt", "w");
    for (long i = 0; i < countnum; i++) {
        fprintf(outres, "%e %e\n", Vrex[i], Irex[i]);
    }
    fclose(outres);
    return;
}

double ChiSq(long countnum, double *Vnum, double *Inum, double *Irex)
{ /*ChiSquare figure of merit calculation*/
    double sum = 0;
    for (long i = 0; i < countnum; i++) {
        sum += (Inum[i] - Irex[i]) * (Inum[i] - Irex[i]) / Irex[i] / Irex[i];
    }
    return sum / countnum;
}

double ChiSqHi(long countnum, double *Vnum, double *Inum, double *Irex)
{ /*ChiSquare figure of merit calculation*/
    double sum = 0;
    for (long i = 0; i < countnum; i++) {
        sum = sum + (Inum[i] - Irex[i]) * (Inum[i] - Irex[i]) * i * 1e8;
    }
    return sum / countnum;
}

double ChiSqDer(long countnum, double *Vnum, double *Inum, double *Irex)
{ /*ChiSquare figure of merit calculation*/
    double sum = 0;
    for (long i = 0; i < countnum - 1; i++) {
        sum += (Inum[i] - Irex[i]) * (Inum[i] - Irex[i]) / Irex[i] / Irex[i]
               + ((Inum[i + 1] - Inum[i]) / (Vnum[i + 1] - Vnum[i])
                  - (Irex[i + 1] - Irex[i]) / (Vnum[i + 1] - Vnum[i]))
                     * ((Inum[i + 1] - Inum[i]) / (Vnum[i + 1] - Vnum[i])
                        - (Irex[i + 1] - Irex[i]) / (Vnum[i + 1] - Vnum[i]))
                     / ((Irex[i + 1] - Irex[i]) / (Vnum[i + 1] - Vnum[i]) * (Irex[i + 1] - Irex[i])
                        / (Vnum[i + 1] - Vnum[i]));
    }
    return sum / countnum;
}

void GetDBStartPoint(double *par, const char *fname)
{ /*Not used for fitting now*/
    CFoo foo(0);
    FILE *DB = fopen(fname, "r");
    int num, numentries = 0, flag = 0, indmin = 0;
    double tmp1, tmp2, tmp3, chimin = 0;
    fscanf(DB, "%d", &num);
    while (flag != EOF) {
        flag = fscanf(DB, "%lf %lf %lf", &tmp1, &tmp2, &tmp3);
        for (int i = 0; i < num; i++)
            flag = fscanf(DB, "%lf %lf", &tmp1, &tmp2);
        if (flag != EOF)
            numentries++;
    }
    freopen(fname, "r", DB);
    fscanf(DB, "%d", &num);
    double *arrI, *arrV, *arrpar1, *arrpar2, *arrpar3;
    arrI = new double[num * numentries];
    arrV = new double[num * numentries];
    arrpar1 = new double[numentries];
    arrpar2 = new double[numentries];
    arrpar3 = new double[numentries];
    for (int i = 0; i < numentries; i++) {
        fscanf(DB, "%lf %lf %lf", &arrpar1[i], &arrpar2[i], &arrpar3[i]);
        for (int j = 0; j < num; j++)
            fscanf(DB, "%lf %lf", &arrV[i * num + j], &arrI[i * num + j]);
    }
    printf("DB loaded successfully: %d entries, %d ppe\n", numentries, num);
    fclose(DB);
    double *tmpI, *tmpV;
    tmpI = new double[num];
    tmpV = new double[num];

    Resample(num, foo.countexp, tmpI, tmpV, foo.Iexp, foo.Vexp, arrI, arrV);
    chimin = ChiSq(num, arrV, arrI, tmpI);
    for (int i = 1; i < numentries; i++) {
        Resample(num, foo.countexp, tmpI, tmpV, foo.Iexp, foo.Vexp, arrI + num * i, arrV + num * i);
        tmp1 = ChiSq(num, arrV + num * i, arrI + num * i, tmpI);
        if (tmp1 < chimin) {
            chimin = tmp1;
            indmin = i;
        }
    }
    par[0] = arrpar1[indmin];
    par[1] = arrpar2[indmin];
    par[2] = arrpar3[indmin];
    printf("Found starting point from DB: %lf %lf %lf, chisq: %e \n", par[0], par[1], par[2], chimin);
    delete[] arrI;
    delete[] arrV;
    delete[] arrpar1;
    delete[] arrpar2;
    delete[] arrpar3;
    delete[] tmpI;
    delete[] tmpV;
    return;
}
