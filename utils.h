#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <string>
#include <valarray>

std::tuple<std::valarray<double>, std::valarray<double>> getExperimentalData(const std::string& filename,
                                                                             bool removeOffset);

std::tuple<std::valarray<double>, std::valarray<double>> Resample(const std::valarray<double>& Iexp,
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

void writeConverg(double x, double f, std::chrono::time_point<std::chrono::steady_clock> start);

#endif
