#ifndef _UTILS_H_
#define _UTILS_H_

long GetExp(char *fname, int iNumCols, int iColNum, double *&Iexp, double *&Vexp, bool bRemOffset);
void Resample(long countnum,
              long countexp,
              double *&Irex,
              double *&Vrex,
              double *&Iexp,
              double *&Vexp,
              double *Inum,
              double *Vnum);
double ChiSq(long countnum, double *Vnum, double *Inum, double *Irex);
double ChiSqHi(long countnum, double *Vnum, double *Inum, double *Irex);
double ChiSqDer(long countnum, double *Vnum, double *Inum, double *Irex);
double IntSq(long countexp, double *Vofx, double *Iofx);
void GetDBStartPoint(double *par, const char *fname);

#endif
