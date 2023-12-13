#ifndef COOF_H
#define COOF_H

#include <valarray>

class OffsetEliminator {
    std::valarray<double> Iref, Vref;

public:
    OffsetEliminator(const std::valarray<double>& Iexp, const std::valarray<double>& Vexp);

    double operator()(double offset) const;
};

#endif
