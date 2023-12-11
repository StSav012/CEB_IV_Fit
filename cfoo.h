#ifndef CFOO_H
#define CFOO_H

#include <valarray>

class CFoo {
    const size_t iNumParams = 21;

    std::valarray<double> Iexp, Vexp;
    std::valarray<double> Inum, Vnum;
    std::valarray<double> Irex, Vrex;

    std::valarray<double> par;
    std::valarray<bool> ToFit;
    size_t iParNum;

public:
    double operator()(double dParam);

    size_t CEB_2eq_parallel_lite();

    void SeqFit(size_t iRunCount);

    explicit CFoo(size_t parnum);
};

#endif
