#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include "nr.h"
#include <string.h>
#include <time.h> 
#include <omp.h>
#include "mins.h"
#include <limits>

long GetExp(char* fname, int iNumCols, int iColNum, double *&Iexp, double *&Vexp, bool bRemOffset);
void Resample(long countnum, long countexp, double *&Irex, double *&Vrex,  double *&Iexp, double *&Vexp, double *Inum, double *Vnum);
double ChiSq(long countnum, double *Vnum, double* Inum, double* Irex);
double ChiSqHi(long countnum, double *Vnum, double* Inum, double* Irex);
double ChiSqDer(long countnum, double *Vnum, double* Inum, double* Irex);
double IntSq(long countexp, double *Vofx, double* Iofx);
void GetDBStartPoint(double* par, const char *fname);

class CFoo {
public:
	long countexp, countnum;
	double *Iexp, *Vexp, *Inum, *Vnum;
	double *Irex, *Vrex;
	double dMinimize;
	double dVStartVg, dVFinVg, dV;
	
	static const int iNumParams=21;
	double par[iNumParams];
	bool ToFit[iNumParams];
	int iParNum;
	struct parametername;
	
	double Pbg; //power for all structure

	double dPbg; //power for 1 bolometer

	double beta; //returned power

	int TephPOW; //returned power

	double gamma; //gap smearing, returned power

	double Vol; //volume of absorber in um^3

	double VolS; //volume of superconductor

	double Z; // heat exchange in normal metal in nW/(K^5*micron^3)

	double ZS; //sigma and volume of supercond-r

	double Tc; //critical temp of superconductor

	double Rn; //normal resistance of 1 SIN, Ohm 

	double Ra; //resistance of absorber, Ohm
	
	int M; //number of bolometers in series,

	int MP; // number of bolometers in parallel

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
class COff {
public: 
	double *Iofx, *Vofx, *Iref, *Vref;
	long count;
	COff(double *Iexp, double *Vexp, long countexp);
	~COff();
	double operator()(double dOffset);
};


