#ifndef _COOF_H_
#define _COOF_H_

class COff
{
public:
    double *Iofx, *Vofx, *Iref, *Vref;
    long count;
    COff(double *Iexp, double *Vexp, long countexp);
    ~COff();
    double operator()(double dOffset);
};

#endif
