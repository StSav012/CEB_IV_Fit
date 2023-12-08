#ifndef _CFOO_H_
#define _CFOO_H_

class CFoo
{
public:
    long countexp, countnum;
    double *Iexp, *Vexp, *Inum, *Vnum;
    double *Irex, *Vrex;
    double dMinimize;
    double dVStartVg, dVFinVg, dV;

    static const int iNumParams = 21;
    double par[iNumParams];
    bool ToFit[iNumParams];
    int iParNum;
    struct parametername;

    double Pbg; //power for all structure

    double dPbg; //power for 1 bolometer

    double beta; //returned power

    int TephPOW; //exponent for Te-ph

    double gamma; //gap smearing, returned power

    double Vol; //volume of absorber

    double VolS; //volume of superconductor

    double Z; //heat exchange in normal metal

    double ZS; //sigma of superconductor

    double Tc; //critical temperature of superconductor

    double Rn; //normal resistance of bolometer

    double Ra; //resistance of absorber

    int M; //number of bolometers in series

    int MP; //number of bolometers in parallel

    double Tp; //phonon temperature

    double Rleak; //leakage resistance of SIN

    float Wt; //transparency of barter

    float tm; //depairing energy

    float ii; //coefficient for Andreev current

    double operator()(double dParam);
    long CEB_2eq_parallel_lite(void);
    void SeqFit(int iRunCount);
    CFoo(int parnum);
    ~CFoo();
};

#endif
