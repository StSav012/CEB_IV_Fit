#ifndef _CFOO_H_
#define _CFOO_H_

#include <valarray>

class CFoo
{
    const size_t iNumParams = 21;

    std::valarray<double> Iexp, Vexp;
    std::valarray<double> Inum, Vnum;
    std::valarray<double> Irex, Vrex;

    std::valarray<double> par;
    std::valarray<bool> ToFit;
    int iParNum;

public:
    double operator()(double dParam);
    long CEB_2eq_parallel_lite(void);
    void SeqFit(int iRunCount);
    CFoo(int parnum);
};

#endif
