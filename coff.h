#ifndef _COOF_H_
#define _COOF_H_

#include <valarray>

class COff
{
    std::valarray<double> Iofx, Vofx, Iref, Vref;

public:
    COff(const std::valarray<double> &Iexp, const std::valarray<double> &Vexp);
    double operator()(double dOffset);
};

#endif
