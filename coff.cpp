#include <cmath>
#include <cstdio>

#include "coff.h"
#include "utils.h"

COff::COff(double *Iexp, double *Vexp, long countexp)
{
    count = countexp;
    Iofx = new double[countexp];
    Vofx = new double[countexp];
    Iref = new double[countexp];
    Vref = new double[countexp];
    for (long i = 0; i < countexp; i++) {
        Vref[i] = Vexp[i];
        Iref[i] = Iexp[i];
    }
}

COff::~COff()
{
    if (Iofx != NULL)
        delete[] Iofx;
    if (Vofx != NULL)
        delete[] Vofx;
    if (Iref != NULL)
        delete[] Iref;
    if (Vref != NULL)
        delete[] Vref;
}

double COff::operator()(double dOffset)
{
    long lLowI = 0;
    double tmp = std::abs(Iref[0]);
    for (long i = 0; i < count; i++) {
        Vofx[i] = Vref[i] - dOffset;
        Iofx[i] = Iref[i];
        if (std::abs(Iref[i]) < tmp) {
            lLowI = i;
            tmp = std::abs(Iref[i]);
        }
    }
    if (Iofx[lLowI] < 0)
        lLowI++;
    long lLengthPos = count - lLowI;
    long lLengthNeg = count - lLengthPos;
    double *Vlow, *Ilow, *Vhigh, *Ihigh;
    long lLengthLow = 0, lLengthHigh = 0;
    if (lLengthPos > lLengthNeg) {
        lLengthHigh = lLengthPos;
        lLengthLow = lLengthNeg;
        Vhigh = new double[lLengthHigh];
        Ihigh = new double[lLengthHigh];
        Vlow = new double[lLengthLow];
        Ilow = new double[lLengthLow];
        for (long i = 0; i < lLengthHigh; i++) {
            Vhigh[i] = std::abs(Vofx[lLowI + i]);
            Ihigh[i] = std::abs(Iofx[lLowI + i]);
        }
        for (long i = 0; i < lLengthLow; i++) {
            Vlow[i] = std::abs(Vofx[lLowI - 1 - i]);
            Ilow[i] = std::abs(Iofx[lLowI - 1 - i]);
        }
    } else {
        lLengthHigh = lLengthNeg;
        lLengthLow = lLengthPos;
        Vhigh = new double[lLengthHigh];
        Ihigh = new double[lLengthHigh];
        Vlow = new double[lLengthLow];
        Ilow = new double[lLengthLow];
        for (long i = 0; i < lLengthHigh; i++) {
            Vhigh[i] = std::abs(Vofx[lLowI - 1 - i]);
            Ihigh[i] = std::abs(Iofx[lLowI - 1 - i]);
        }
        for (long i = 0; i < lLengthLow; i++) {
            Vlow[i] = std::abs(Vofx[lLowI + i]);
            Ilow[i] = std::abs(Iofx[lLowI + i]);
        }
    }
    double *Irex = new double[lLengthLow], *Vrex = new double[lLengthLow];
    FILE *test1 = fopen("test1.txt", "w+");
    for (int i = 0; i < lLengthLow; i++)
        fprintf(test1, "%e %e \n", Vlow[i], Ilow[i]);
    fclose(test1);
    FILE *test2 = fopen("test2.txt", "w+");
    for (int i = 0; i < lLengthHigh; i++)
        fprintf(test2, "%e %e \n", Vhigh[i], Ihigh[i]);
    fclose(test2);
    Resample(lLengthLow, lLengthHigh, Vrex, Irex, Vhigh, Ihigh, Vlow, Ilow);
    FILE *test3 = fopen("test3.txt", "w+");
    for (int i = 0; i < lLengthLow; i++)
        fprintf(test3, "%e %e \n", Vrex[i], Irex[i]);
    fclose(test3);
    double dResult = ChiSqHi(lLengthLow, Ilow, Vlow, Vrex);
    if (Ilow != NULL)
        delete[] Ilow;
    if (Vlow != NULL)
        delete[] Vlow;
    if (Ihigh != NULL)
        delete[] Ihigh;
    if (Vhigh != NULL)
        delete[] Vhigh;
    if (Irex != NULL)
        delete[] Irex;
    if (Vrex != NULL)
        delete[] Vrex;
    return dResult;
}
