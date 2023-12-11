#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <valarray>

size_t GetExp(const std::string& fname,
              std::valarray<double>& Iexp,
              std::valarray<double>& Vexp,
              bool bRemOffset);

void Resample(std::valarray<double>& Irex,
              std::valarray<double>& Vrex,
              const std::valarray<double>& Iexp,
              const std::valarray<double>& Vexp,
              const std::valarray<double>& Inum,
              const std::valarray<double>& Vnum);

double ChiSq(const std::valarray<double>& Inum, const std::valarray<double>& Irex);

double ChiSqHi(const std::valarray<double>& Inum, const std::valarray<double>& Irex);

double ChiSqDer(const std::valarray<double>& Vnum,
                const std::valarray<double>& Inum,
                const std::valarray<double>& Irex);

// void GetDBStartPoint(double *par, const char *fname);

void writeIV(const std::string& filename,
             const std::valarray<double>& I,
             const std::valarray<double>& V);

void writeConverg(double x, double f, time_t start);

#endif
