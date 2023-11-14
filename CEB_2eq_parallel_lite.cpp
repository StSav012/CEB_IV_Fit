// CEB_2eq.cpp : Defines the entry point for the console application.
//Pleak is included in Pcool

//for bolometers (absorber + 2SINs) w/ Andreev Current

//#include <complex.h>
#include "header.h"

#define PI 3.14159265358979
#define ME 9.10938188e-31	//kg
#define EL 1.60217646e-19	//coulombs
#define HPLANCK 6.626068e-34	//m2 kg / s
#define HBAR 1.05457148e-34	//m2 kg / s
#define KBOL 8.6173324e-5	//eV K-1
#define THREADS 56

//1 electron volt = 1.60217646e-19 joule

double current(float v, float tau)
{//a procedure which computes the current using approximation
	//i is measured in [delta(0) / eRn]
	/*'v' is voltage, 'tau' is for temperature*/
	
	double i;
	double a0, a1, a2, a3;
	float s;

	a0 = 1 + 0.375 * tau - 0.1171875 * tau * tau;
	a1 = tau * a0 * a0;
	a2 = 1 + exp((abs(v) - 1) / tau - (1.15 + tau));
	a3 = (v * v - 1) / (1 + exp(-(v * v - 1) / tau));

	if (v > 0)
		s = 1;
	else s = -1;

    	i = s * sqrt(2 * PI * a1 / a2 + a3) * (1 / (2 * exp((1 - v) / tau) + 1) + 1 / (2 * exp((1 + v) / tau) + 1));
	return i;
}

double currentInt(float DT, float v, float tau, float tauE)
{//computes the current using exact integral 
	/*'DT' is delta (energy gap), 'v' is voltage, 'tau' is for temperature, 'tauE' is for electron temperature*/
	double i, v1;
	double a0, a1, a2; 
	double x; //energy
	float dx;
	int l, inf;
	float s;
	dx = 0.2e-3;
	inf = (int)(20 / dx);
	i = 0;
	x = DT + dx;
	v1 = abs(v);

	double prev_i = 0;     // Previous 'i' value
	double accuracy = 1e-6;  // Calculation accuracy

//#pragma omp parallel for default(shared) private(x, a0, a1, a2) reduction(+:i) num_threads(THREADS)
	for (l = 0; l <= inf; l++)
	{  
		a0 = 1 / (exp((x - v1) / tauE) + 1);
		a1 = 1 / (exp(x / tau) + 1); 
	 	a2 = sqrt(x * x - DT * DT);
		i += abs(x) / a2 * (a0 - a1);
		x += dx;

		// Accuracy check on every iteration
		if (fabs(i - prev_i) < accuracy)
			break;
		prev_i = i;
	}
	x= -DT - dx;

//#pragma omp parallel for default(shared) private(x, a0, a1, a2) reduction(+:i) num_threads(THREADS)
   	for (l = 0; l <= inf; l++)
   	{  
		a0 = 1 / (exp((x - v1) / tauE) + 1);
		a1 = 1 / (exp(x / tau) + 1);
		a2 = sqrt(x * x - DT * DT);
		i += abs(x) / a2 * (a0 - a1);
		x -= dx;
	   
		// Accuracy check on every iteration
		if (fabs(i - prev_i) < accuracy)
			break;
		prev_i = i;
	}
	if (v > 0)
		s = 1;
	else s = -1;
	i = i * dx * s;
	return i;
}

double PowerCool(float v, float tau, float tauE)
{//a procedure which computes the cooling power of SIN junction 
   
	double p;
	double a0, a1, a2, a3;
	double b0, b1, b2, b3;
	double c0, c1;

	a3 = 0; 
	b3 = 0;
	v = abs(v);

	a0 = (1 - v) / tauE;
	a1 = sqrt(2 * PI * tauE) * ((1 - v) / (2 * exp(a0) + 1.28) + 0.5 * tauE / (2 * exp(a0) + 0.64)) / (exp(-2.5 * (a0 + 2)) + 1);
	if ((v - 1 - tauE) > 0) 
	{
		a2 = sqrt(v * v - 1);
		a3 = 0.5 * (-v * a2 + log(v + a2) + PI * PI * tauE * tauE / 3 * v / a2) / (exp(2.5 * (a0 + 2)) + 1);
	}

	b0 = (1 + v) / tauE;
	b1 = sqrt(2 * PI * tauE) * ((1 + v) / (2 * exp(b0) + 1.28) + 0.5 * tauE / (2 * exp(b0) + 0.64)) / (exp(-2.5 * (b0 + 2)) + 1);
	if ((-v - 1 - tauE) > 0) 
	{
		b2 = sqrt(v * v - 1);
		b3 = 0.5 * (v * b2 + log(-v + b2) - PI * PI * tauE * tauE / 3 * v / b2) / (exp(2.5 * (b0 + 2)) + 1);
	}

	c0 = 1 / tau;
	c1 = 2 * sqrt(2 * PI * tau) * (1 / (2 * exp(c0) + 1.28) + 0.5 * tau / (2 * exp(c0) + 0.64)) / (exp(-2.5 * (c0 + 2)) + 1);
   
	p = (a1 + a3 + b1 + b3 - c1);
	return p;
}

double AndCurrent(float DT, float v, float tauE, float Wt, float tm)
{//computes the Andreev current using exact integral 
	//Wt - omega with ~ from Vasenko 2010
	//dd - 1/tm, energy of pair breaking
   
	double i, v1;
	double a0, a1, a2; 
	double x; //energy
	float dx;
	int l, inf;
	//FILE *f9;
  
	dx = 0.2e-3;
	inf = (int)(DT / dx) - 1;

	i = 0;
	x = dx;
	v1 = abs(v);

	//f9=fopen("SINi.dat","w");
	
//#pragma omp parallel for default(shared) private(a0, a1, a2) reduction(+:i) num_threads(THREADS)
	for (l = 0; l <=inf-1; l++)
	{  
		a0 = tanh((x + v) / (2 * tauE));
		a1 = tanh((x - v) / (2 * tauE));
		a2 = 2 * Wt * sqrt(DT * DT - x * x) / tm / (pow(2 * Wt * x - x * sqrt(DT * DT - x * x) / DT, 2) + (DT * DT - x * x) / pow(DT * tm, 2));
		i += DT / sqrt(DT * DT - x * x) * (a0 - a1) * a2;
		//fprintf(f9,"%g %g\n", x, i); 
		x += dx;
	}
   	i = i * dx;
	return i;
	//fclose(f9); 
}

double PowerCoolInt(float DT, float v, float tau, float tauE, double *Ps)
{//integral of the cooling power of SIN junction 
	//E = x, V = v, T = tau, [x] = [v] = [tau]
   
	double p, q;
	double a, a1, a2;
	double b, b1;
	double x, dx, de;
	int l, inf;
	double po;	//power   
	FILE *f8;	//density of states
   
	po = 0;
	*Ps = 0;
	v = abs(v);

	dx = 2e-4;
	inf = (int)(20 / dx);

	//f8 = fopen("DOS.txt","w");
  
	de = DT;
	x = de + dx;

//#pragma omp parallel for default(shared) private(x, a, a1, a2) reduction(+:po) num_threads(THREADS)
	for (l = 0; l <= inf; l++)
	{  
		a = sqrt(x * x - DT * DT);	   

		a2 = 1 / (exp((x - v) / tauE) + 1);
		a1 = 1 / (exp(x / tau) + 1);

		po += abs(x) * (x - v) * (a2 - a1) / a;	//po += abs(x) * (x - v) / a * (a2 - a1);
		*Ps += abs(x) * (x) * (a2 - a1) / a;
		x += dx;
	}
	x = -de - dx;
  
//#pragma omp parallel for default(shared) private(x, a, a1, a2) reduction(+:po) num_threads(THREADS)
	for (l = 0; l <= inf; l++)
	{  
		a = sqrt(x * x - DT * DT);	   

		a2 = 1 / (exp((x - v) / tauE) + 1);
		a1 = 1 / (exp(x / tau) + 1);

		po += abs(x) * (x - v) * (a2 - a1) / a;
		*Ps += abs(x) * x * (a2 - a1) / a;
		x -= dx;
	}
	po = po * dx;
	*Ps = *Ps * dx;
	return po;
}

long CFoo::CEB_2eq_parallel_lite(void)
{
	FILE *f2, *f3, *f3old, *f4, *f5;					//IV 
	
	int k,q,j,l,n;
	
	int Nt, NV;

	double I[10002], I1[10002];
	
	float I_A[10002], I_As[10002];
	
	double V[10002];
	
	double Ib;								//bias for electron cooling

	double I0;								//units of current

	double Vstr, Vfin;							//range of voltage

	float Z;								//heat exchange in normal metal in nW / (K^5 * micron^3)

	double Delta;

	double E, dE, E1;							//energy, eV

	double int1, int2;							//integrals

	double Pbg;								//background power, pW

	double Te, dT, Ts;							//electron and phonon temperature

	double tau, Vg, tauE, tauC, T1, T2;					//dimentionless variables
	
	double tauold;

	double Pe_p, Pabs, Pleak, Pheat, Pcool, P, Pcool1, Ps, Ps1, Pand;	//for power

	double Rsg, Rsg1;							//subgap resistance of 1 SIN, kOhm 

	double Rleak;								//leakage resistance of SIN

	double Rsin;								//normal resistance of single SIN junction

	double Rn1, Rleak1, V1, V2, Vcur;					//for thermometer junctions

	double G, dPdT, dIdT, dIdV, dPdV, Sv, dPT, mm;				//for noise 
	
	double dP, dI, dPdI;
	
	double NEP, NEPs, NEPep, NEPa, NEPph;
	
	double vn, in;								//amplifier voltage and current noise

	double Tp0, x, DeltaT;
	
	double eps;								//coefficient between Ps and Ts

	double NoiA;								//amplifier noise

	float G_NIS, G_e; 
	
	float Wt;								//transparency of barrier

	float tm;								//depairing energy

	double I_A0, I_As0;

	float ii;								//coefficient for Andreev current

	clock_t start, finish;
	
	start = clock();
	
	char c = 0;
      
// --------- known/guessed physical parameters

	Pbg = par[0]; 			//incoming power for all structure, pW

	beta = par[1];			//returning power ratio, <1

	TephPOW = par[2];		//exponent for Te-ph

	gamma = par[3];			//gap smearing, not used, should be 0

	Vol = par[4];			//volume of absorber, um^3

	VolS = par[5];			//volume of superconductor, um^3

	Z = par[6];			//heat exchange in normal metal, nW/(K^5*um^3)

	ZS = par[7];			//sigma

	Tc = par[8];			//critical temperature, K

	Rn = par[9] / M * MP;		//normal resistance for 1 bolometer, Ohm

	Rleak = par[10] / M * MP;	//leakage resistance for 1 bolometer, Ohm

	Wt = par[11];			//transparency of barrier

	tm = par[12];			//depairing energy

	ii = par[13];			//coefficient for Andreev current

	Ra = par[14];			//normal resistance for 1 absorber, Ohm

	M = par[15];			//number of bolometers in series

	MP = par[16];			//number of bolometers in parallel

	Tp = par[17];			//phonon temperature, K

	dVFinVg = par[18];		//voltage range end

	dVStartVg = par[19];		//voltage range start

	dV = par[20];			//voltage step, V

	Te = Tp;			//electron temperature, is to be found, K

	Ts = Tp;			//electron temperature in superconductor, K
	
	//Ib = 1.4;			//bias current for electron cooling, nA

	dPbg = Pbg / M / MP;		//incoming power for 1 bolometer, pW

	Delta = 1.764 * Tc;		//in K, Vg[eV] = Tc * 1.764 * 86.25e-6
/*
 	char* cPbgNoise = new char[6];
	_gcvt_s(cPbgNoise, sizeof(cPbgNoise), Pbg, 5);
	strcat(cPbgNoise, " pW Noise.dat");

	char* cPbgTe = new char[6];
	_gcvt_s(cPbgTe, sizeof(cPbgTe), Pbg, 5);
	strcat(cPbgTe, " pW Te.txt");

	char* cPbgNEP = new char[6];
	_gcvt_s(cPbgNEP, sizeof(cPbgNEP), Pbg, 5);
	strcat(cPbgNEP, " pW NEP.dat");

	char* cPbgG = new char[6];
	_gcvt_s(cPbgG, sizeof(cPbgG), Pbg, 5);
	strcat(cPbgG, " pW G.txt");
 */
	f3old = fopen("Te_old.txt", "w");
	f3 = fopen("Te.txt", "r");
	//f3 = fopen(cPbgTe, "r");
	
	if (f3 != NULL)
	{
		while (fscanf(f3, "%c", &c) != EOF) 
			fprintf(f3old, "%c", c);
		fclose(f3);
	}
	
	fclose(f3old);
/*
 	f2 = fopen(cPbgNoise, "w");
	f3 = fopen(cPbgTe, "w"); 
	f4 = fopen(cPbgNEP, "w");
	f5 = fopen(cPbgG, "w");
 */
	f2 = fopen("Noise.dat", "w");
	f3 = fopen("Te.txt", "w"); 
	f4 = fopen("NEP.dat", "w");
	f5 = fopen("G.txt", "w");

	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);

//---------- normalized constants

	Rn = (Rn - Ra) / 2;
	//Rleak=Rn/gamma;	//Ohm, is not needed if Gamma<>0

	I0 = Delta / Rn * KBOL * 1e9;
	
	Vg = Delta / 11604.505;
	/*Vg = Delta * KBOL;*/  
	
	printf("%lf\n", Vg);

	tau = Ts / Delta;
	
	tauE = Te / Delta;
	
	tauC = Tc / Delta;

//---------- calculation parameters

	//voltage step, V
	
	Vfin = dVFinVg * Vg;
	
	Vstr = dVStartVg * Vg;
	
	NV = (Vfin - 0.0 * Vg) / dV;
	
	if (Inum == NULL) 
		Inum = new double[NV - 1];
	
	if (Vnum == NULL)  
		Vnum = new double[NV - 1];
	
	for (j = 0; j <= NV; j++) 
		V[j] = Vstr + ((j - 0) * dV);	// = [V]

	Tp0 = Tp;
	
	DeltaT = 1;
/*
	f2 = fopen(cPbgNoise, "a");
	f3 = fopen(cPbgTe, "a"); 
  	f4 = fopen(cPbgNEP, "a");
	f5 = fopen(cPbgG, "a");
*/
	f2 = fopen("Noise.dat", "a");
	f3 = fopen("Te.txt", "a"); 
  	f4 = fopen("NEP.dat", "a");
	f5 = fopen("G.txt", "a");

	fprintf(f2, "Voltage\tNOISEep\tNOISEs\tNOISEa\tNOISE\tNOISEph\tNOISE^2-NOISEph^2\n");
	fprintf(f3, "Voltage\tCurrent\tIqp\tIand\tV/Rleak\tTe\tTs\tDeltaT\tPand\tPleak\tPabs\tPcool\n");
	fprintf(f4, "Voltage\tCurrent\tNEPep\tNEPs\tNEPa\tNEP\tNEPph\tSv\tNEP^2-NEPph^2\n");
	fprintf(f5, "Voltage\tGe\tGnis\n");

	FILE* check=fopen("conv.txt", "w");
	
	for (j = 1; j <= NV - 1; j++)   //next voltage
	{
		for (n = 1; n <= 5; n++)   //next interation
		{			
			T1 = 0; 
			
			T2 = 3 / 1.764;
			
			tauE = (T1 + T2) / 2;
			
			DeltaT = sqrt(1 - pow(Ts / Tc, float(3.2))); 
			
			for (l = 1; l <= 15; l++) //next iteration
			{
				tauE = (T1 + T2) / 2;
					
				I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9;	//in nA
					
				I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0;	//in nA
					
				Pe_p = Z * Vol * (pow(Tp, TephPOW/*7,6,5*/) - pow(tauE * Delta, TephPOW/*7,6,5*/)) * 1e3;	//pW
					
				Pabs = I[j] * I[j] * Ra * 1e-6;	//pW
					
				Pleak = 2.0 * V[j] * V[j] / Rleak * 1e12;	// pW, 2 because of 2 SIN
					
				Pand = (I_A[j] * 1e-3) /*uA*/ * (I_A[j] * 1e-3) /*uA*/ * Ra /*Ohm*/ + 2.0 * (I_A[j] * 1e3) /*pA*/ * V[j] /*V*/;	//pW, absorber + Andreev
					
				Pcool = 1.0 * PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE, &Ps) * Vg * Vg / Rn * 1e12;	//pW, 0.1 - experimental parameter
					
				Ps = Ps * Vg * Vg / Rn * 1e12;	//returning power from S to N
					
				Pheat = Pe_p + 1.0 * Pabs + 1.0 * Pand + dPbg + 2.0 * beta * Ps + 1.0 * Pleak;	// + Pand;
					
				P = Pheat - 2.0 * Pcool;
					
				if (P < 0)
				{
					T2 = tauE;
				}	
						
				else 
				{
					T1 = tauE;
				}
			}
		}
		//DeltaT = 1; // !!!
		
		Te = tauE * Delta;
		
		Ts = tau * Delta;
		
		DeltaT = sqrt(1 - pow(Ts / Tc, float(3.2)));

		I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9;
		
		I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0;
	   
		fprintf(f3, "%f\t%g\t%g\t%g\t%f\t%f\t%g\t%g\t%g\t%g\t%g\t%g\n", M * (2 * V[j] + (I[j] * 1e-9 + 1.0e-9 * I_A[j]) * Ra), MP * 1e-9 * (I[j] + I_A[j]), I[j] * MP * 1e-9, I_A[j] * MP * 1e-9, (V[j] / Rleak * 1e9) * MP, Te, Ts, DeltaT, Pand, Pleak, Pabs, Pcool);
		
		printf("Voltage: %g Current: %g\n", M * (2 * V[j] + (I[j] * 1e-9 + 1.0e-9 * I_A[j]) * Ra), MP * 1e-9 * (I[j] + I_A[j]) /*, I[j] * MP * 1e-9, I_A[j] * MP * 1e-9, Pand, Pleak*/ ); 

		//fprintf(f4,"%f %f\n", M * (2 * V[j] + I[j] * Ra * 1e-9), 2 * Rsg);  
		
		Inum[j - 1] = MP * 1e-9 * (I[j] + I_A[j]);
		
		Vnum[j - 1] = M * (2 * V[j] + ((I[j] + I_A[j]) * Ra * 1e-9));

//----- NEP -------------------------------------
//308mK_1.txt SampleC_200mK_48bolo.txt OL65_305mK_0.txt OL65 OPA111.txt OL76_304.txt
		
		dT = 0.005;

		vn = (3.2e-9) * sqrt(2.0);	//V/sqrt(Hz)	//for 2 amp * sqrt(2) AD745
		in = (6.9e-15) * 1 / sqrt(2.0);	//A/sqrt(Hz)	//for 2 amp * 1 / sqrt(2)

//		vn = (8.0e-9) * sqrt(2.0);	//V/sqrt(Hz)	//for 2 amp * sqrt(2) OPA111
//		in = (0.8e-15) * 1 / sqrt(2.0);	//A/sqrt(Hz)	//for 2 amp * 1 / sqrt(2)

//		vn = (0.9e-9) * sqrt(2.0);	//V/sqrt(Hz)	//for 2 amp * sqrt(2) AD797
//		in = (2.0e-12) * 1 / sqrt(2.0);	//A/sqrt(Hz)	//for 2 amp * 1 / sqrt(2)

//		vn = (1.1e-9) * sqrt(2.0);	//V/sqrt(Hz)	//for 2 amp * sqrt(2) IFN146
//		in = (0.3e-15) * 1 / sqrt(2.0);	//A/sqrt(Hz)	//for 2 amp * 1 / sqrt(2)

//		vn = (5.1e-9) * sqrt(2.0);	//V/sqrt(Hz)	//for 2 amp * sqrt(2) OPA1641
//		in = (0.8e-15) * 1 / sqrt(2.0);	//A/sqrt(Hz)	//for 2 amp * 1 / sqrt(2)

		dPT = PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta, &Ps) - PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta, &Ps);
		
        	dPdT = Vg * Vg / Rn * 1e12 * (dPT) / (2 * dT);	//pW/K
		
	    	dIdT = I0 * (currentInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta) - currentInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta)) / (2 * dT);	//nA/K
		
	    	dIdV = (I0 * (currentInt(DeltaT, V[j + 1] / Vg, tau, tauE) + ii * AndCurrent(DeltaT, V[j + 1] / Vg, tauE, Wt, tm) - currentInt(DeltaT, V[j - 1] / Vg, tau, tauE) - ii * AndCurrent(DeltaT, V[j - 1] / Vg, tauE, Wt, tm))) / (2 * dV);	//nA/V
		
	    	dPdV = Vg * Vg / Rn * 1e12 * (PowerCoolInt(DeltaT, V[j + 1] / Vg, tau, tauE, &Ps) - PowerCoolInt(DeltaT, V[j - 1] / Vg, tau, tauE, &Ps)) / (2 * dV);	//pW/V

        	G_NIS = dPdT; G_e = 5 * Z * Vol * pow(Te, 4) * 1e3;	//pW/K

		G = G_e + 2 * (G_NIS - dIdT / dIdV * dPdV);	//pW/K (2)SINs		
	   
	    	Sv = - 2 * dIdT / dIdV / G;	//V/pW  for 1 bolo

		Sv = Sv / MP;
		
	    	NEPep = 10 * EL * KBOL * Z * Vol * (pow(Tp, TephPOW) + pow(Te, TephPOW)) * 1e3 * 1e12;	//^2 pW^2/Hz

		NoiA = (vn * vn + pow(float(in * (2 * 1e9 / dIdV + Ra) * M / MP), float(2)));	//(V/sqrt(Hz))^2
		
		NEPa = (NoiA) / Sv / Sv;	//pW^2/Hz

//----- NEP SIN approximation ----------------------------

		dI = 2 * EL * abs(I[j]) / pow(dIdV * Sv, 2) * 1e9;	//pW^2/Hz
		
	    	dPdI = 2 * 2 * EL * Pcool / (dIdV * Sv) * 1e9;	//pW^2/Hz, second '2' is from comparison with integral
		
		mm = log(sqrt(2 * PI * KBOL * Te * Vg) / (2 * abs(I[j]) * Rn * 1e-9));
		
	    	dP = (0.5 + mm * mm) * (KBOL * Te) * (KBOL * Te) * abs(I[j]) * EL * 1e-9 * 1e24;	//pW^2/Hz
		
	    	NEPs = 2 * (dI - 2 * dPdI + dP);	//pW^2/Hz (2) //all terms positive
		
		//fprintf(f5, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP); 

//----- NEP SIN integral ----------------------------
/* 
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / pow(dIdV * Sv, 2) * 1e9;	//pW^2/Hz
  
		dPdI = dPdI / (dIdV * Sv) * EL * Vg * Vg / Rn * 1e12 * 1e9;	//pW^2/Hz
  
		dP = dP * Vg * Vg * Vg / Rn * EL * 1e24;	//pW^2/Hz 
        
		NEPs = 2 * (dI - 2 * dPdI + dP);	//pW^2/Hz (2)  
  
		fprintf(f9, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);
*/
//-----------------------------------------------
//-----------------------------------------------
		
		NEPph = sqrt(M * MP * 2 * Pbg * 0 * 1e9 * HPLANCK * 1e-12 + pow(M * MP * Pbg * 1e-12, 2) / 1e3 / 1e9);	//W/sqrt(Hz) at 0 GHz
		
		//NEPph = sqrt(M * MP * 2 * Pbg * 350 * 1e9 * HPLANCK * 1e-12 + pow(M * MP * Pbg * 1e-12, 2) / 1.552 / 1e9);	//W/sqrt(Hz) at 350 GHz
		
		NEPph = NEPph * 1e12;	//pW

		NEP = sqrt(M * MP * (NEPep + NEPs) + NEPa + NEPph * NEPph);	//all squares 

		fprintf(f2, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", M * (2 * V[j] + I[j] * Ra * 1e-9), sqrt(M * MP * NEPep) * 1e9 * abs(Sv), sqrt(M * MP * NEPs) * 1e9 * abs(Sv), sqrt(NoiA) * 1e9, NEP * 1e9 * abs(Sv), NEPph * 1e9 * abs(Sv), 1e9 * abs(Sv) * sqrt(NEP * NEP - NEPph * NEPph)); 

		fprintf(f4, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", M * (2 * V[j] + I[j] * Ra * 1e-9), MP * 1e-9 * (I[j]), sqrt(M * MP * NEPep) * 1e-12, sqrt(M * MP * NEPs) * 1e-12, sqrt(NEPa) * 1e-12, NEP * 1e-12, NEPph * 1e-12, abs(Sv) * 1e12, 1e-12 * sqrt(NEP * NEP - NEPph * NEPph)); 
		
		fprintf(f5, "%g\t%g\t%g\n", M *(2 * V[j] + I[j] * Ra * 1e-9), G_e, G_NIS); 
		
		printf("Sv: %g Te: %g NEPs: %g NEPt: %g\n\n", abs(Sv) * 1e12, Te, sqrt(M * MP * NEPs) * 1e-12, NEP * 1e-12);

		//system("cls"); 
	}

	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);

	finish = clock();
	double sTime = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("\nTime spent: %f sec.", sTime);

	return NV - 1;
}

#undef PI
#undef ME  
#undef E  
#undef HPLANCK
#undef HBAR
#undef KBOL
