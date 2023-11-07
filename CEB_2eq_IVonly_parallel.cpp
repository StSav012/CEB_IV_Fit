// CEB_2eq.cpp : Defines the entry point for the console application.
//Pleak is included in Pcool

//for bolometers (absorber + 2SINs)

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <math.h>
#include <time.h> 
#include <omp.h>
//#include <complex.h>

#define PI 3.14159265358979
#define ME  9.10938188e-31 //kg
#define EL  1.60217646e-19 //coulombs
#define HPLANCK 6.626068e-34 //m2 kg / s
#define HBAR 1.05457148e-34 //m2 kg / s
#define KBOL 8.6173324e-5 //eV K-1
#define THREADS 6

//1 electron volt = 1.60217646 * 10-19 joule



double current(float v, float tau)
{//a procedure which computes the current using approximation
	//i is measured in [delta(0)/eRn]
	double i;
	double a0, a1, a2, a3;
	float s;

	a0 = 1 + 0.375*tau - 0.1171875*tau*tau;
	a1 = tau*a0*a0;
	a2 = 1 + exp((abs(v) - 1)/tau - (1.15+tau));
	a3 = (v*v-1)/(1 + exp(-(v*v-1)/tau));

	if (v>0) s=1;
	else s=-1;

    i = s*sqrt(2*PI*a1/a2 + a3)*(1/(2*exp((1-v)/tau)+1) + 1/(2*exp((1+v)/tau)+1));
	return i;
}

double currentInt(float v, float tau, float tauE)
{//a procedure which computes the current using exact integral 
   
   double i, v1;
   double a0, a1, a2; 
   double x; //energy
   float dx;
   int l, inf;
   float s;

   dx=0.2e-3;
   inf=(int)(20/dx);

   i=0;
   x=1;
   v1=abs(v);
   
   for (l=0; l<=inf; l++)
   {  
	   a0=1/(tauE*(exp((x-v1)/tauE)+ 2 + exp(-(x-v1)/tauE)));
       a1=1/(tau*(exp(x/tau)+ 2 + exp(-x/tau)));
       a2=sqrt(x*x-1);
	   i+=a2*(a0-a1);
	   x+=dx;
   }

   x=-1;
   for (l=0; l<=inf; l++)
   {  
	   a0=1/(tauE*(exp((x-v1)/tauE)+ 2 + exp(-(x-v1)/tauE)));
       a1=1/(tau*(exp(x/tau)+ 2 + exp(-x/tau)));
       a2=sqrt(x*x-1);
	   i-=a2*(a0-a1);
	   x-=dx;
   }

   if (v>0) s=1;
	else s=-1;
   i=i*dx*s;

   return i;
}

double currentIntCom(float DT, float v, float tau, float tauE, float gamma)
{//integral of the current through SIN junction with gap smearing
	//E = x, V=v, T=tau, [x]=[v]=[tau]
	//DT = delta(tau)/delta(0)
   
   double p, q;
   double a, a1, a2;
   double b, b1;
   double x, dx;
   int l, inf;
   double i; //current
   FILE *f8; //density of states
   double tmp1;
   
   i=0;
   //v=abs(v);

   dx=0.2e-3;
   //dx<1e-3 for gamma=1e-4
   //dx<0.2e-3 for gamma=1e-5
   inf=(int)(20/dx);

   //f8=fopen("DOS.txt","w");
  
   x=dx;
#pragma omp parallel for default(shared) private(x, a, b, tmp1, p, q, b1, a2) reduction(+:i) num_threads(THREADS)
	   for (l=0; l<=inf; l++)
	   {  
		   x=dx*(l+1);
		   
		   a=x*x-gamma*gamma-DT;	   b=-2*x*gamma;

		   tmp1=sqrt(a*a+b*b);

		   p=sqrt((a+tmp1)/2); q=sqrt((tmp1-a)/2);

		   b1=(x*p+gamma*q)/(p*p + q*q)*1e6; //fprintf(f8,"%f %f\n",x,b1);

		   a2=1e6/(exp((x-v)/tauE)+ 1); a1=1e6/(exp(x/tau)+ 1);

		   i+=b1*(a2-a1); //fprintf(f8,"%f %f\n",x,b1*(a2-a1));
	   }

   //fprintf(f8,"\n");
   x=-dx;
#pragma omp parallel for default(shared) private(x, a, b, tmp1, p, q, b1, a2) reduction(+:i) num_threads(THREADS)
	   for (l=0; l<=inf; l++)
	   {  
		   x=-dx*(l+1);
		   a=x*x-gamma*gamma-DT;	   b=-2*x*gamma;

		   tmp1=sqrt(a*a+b*b);

		   p=sqrt((a+tmp1)/2); q=-sqrt((tmp1-a)/2);

		   b1=-(x*p+gamma*q)/(p*p + q*q)*1e6; //fprintf(f8,"%f %f\n",x,b1);

		   a2=1e6/(exp((x-v)/tauE)+ 1); a1=1e6/(exp(x/tau)+ 1);

		   i+=b1*(a2-a1); //fprintf(f8,"%f %f\n",x,b1*(a2-a1));
	   }
   //fprintf(f8,"%f\n",i*1e10);
   //fclose(f8);

   //if (v>0) i=i*dx;
   //else i=i*dx*(-1);
   i=i*dx;  //do not forget to add *1e-12!!!
   
   return i;

}

double PowerCool(float v, float tau, float tauE)
{//a procedure which computes the cooling power of SIN junction 
   
   double p;
   double a0, a1, a2, a3;
   double b0, b1, b2, b3;
   double c0, c1;

   a3=0; b3=0;
   v=abs(v);

   a0 = (1-v)/tauE;
   a1 = sqrt(2*PI*tauE)*((1-v)/(2*exp(a0)+1.28)+0.5*tauE/(2*exp(a0)+0.64))/(exp(-2.5*(a0+2))+1);
   if ((v-1-tauE)>0) 
   {
      a2=sqrt(v*v-1);
	  a3=0.5*(-v*a2 + log(v+a2) + PI*PI*tauE*tauE/3*v/a2)/(exp(2.5*(a0+2))+1);
   }

   b0 = (1+v)/tauE;
   b1 = sqrt(2*PI*tauE)*((1+v)/(2*exp(b0)+1.28)+0.5*tauE/(2*exp(b0)+0.64))/(exp(-2.5*(b0+2))+1);
   if ((-v-1-tauE)>0) 
   {
      b2=sqrt(v*v-1);
	  b3=0.5*(v*b2 + log(-v+b2) - PI*PI*tauE*tauE/3*v/b2)/(exp(2.5*(b0+2))+1);
   }

   c0 = 1/tau;
   c1 = 2*sqrt(2*PI*tau)*(1/(2*exp(c0)+1.28)+0.5*tau/(2*exp(c0)+0.64))/(exp(-2.5*(c0+2))+1);
   
   p = (a1 + a3 + b1 + b3 - c1);
   
   return p;

}


double PowerCoolIntCom(float DT, float v, float tau, float tauE, float gamma, double *Ps)
{//integral of the cooling power of SIN junction 
	//E = x, V=v, T=tau, [x]=[v]=[tau]
   
   double p, q;
   double a, a1, a2;
   double b, b1;
   double x, dx;
   int l, inf;
   double po; //power
   double tmPs=0;
   FILE *f8; //density of states
  // double tmp1, tmp2;
   
   po=0;
   *Ps=0;
   v=abs(v);

   dx=0.2e-3;
   inf=(int)(20/dx);

   //f8=fopen("DOS.txt","w");
  
   x=dx;

#pragma omp parallel for default(shared) private(x, a, b, p, q, b1, a2) reduction(+:po, tmPs) num_threads(THREADS)
	   for (l=0; l<=inf; l++)
	   {  
		   x=dx*(l+1);
		   a=x*x-gamma*gamma-DT;	   b=-2*x*gamma;

		   //tmp1=sqrt(a*a+b*b);
		   //p=sqrt((a+tmp1)/2); q=sqrt((tmp1-a)/2);
			p=sqrt((a+sqrt(a*a+b*b))/2); q=sqrt((sqrt(a*a+b*b)-a)/2);

		   b1=(x*p+gamma*q)/(p*p + q*q); //fprintf(f8,"%f %f\n",x,b1);

		   a2=1/(exp((x-v)/tauE)+ 1); a1=1/(exp(x/tau)+ 1);

		  // tmp2=b1*(a2-a1);

		   po+=(x-v)*b1*(a2-a1);
		   tmPs+=x*b1*(a2-a1);
	   }
 
   //fprintf(f8,"\n");
   x=-dx;
   #pragma omp parallel for default(shared) private(x, a, b,  p, q, b1, a2) reduction(+:po, tmPs) num_threads(THREADS)
	   for (l=0; l<=inf; l++)
	   {  
		   x=-dx*(l+1);
		   a=x*x-gamma*gamma-DT;	   b=-2*x*gamma;
		   
		   //tmp1=sqrt(a*a+b*b);
		   //p=sqrt((a+tmp1)/2); q=-sqrt((tmp1-a)/2);
			p=sqrt((a+sqrt(a*a+b*b))/2); q=-sqrt((sqrt(a*a+b*b)-a)/2);
		   b1=-(x*p+gamma*q)/(p*p + q*q); //fprintf(f8,"%f %f\n",x,b1);

		   a2=1/(exp((x-v)/tauE)+ 1); a1=1/(exp(x/tau)+ 1);
		
		   //tmp2=b1*(a2-a1);

		   po+=(x-v)*b1*(a2-a1);
		   *Ps+=x*b1*(a2-a1);
		}
   //fclose(f8);
   
   po=po*dx;
   *Ps=tmPs*dx;
   
   return po;

}

double NEPIntCom(float v, float tau, float tauE, float gamma, double *dI, double *dP, double *dPdI)
{//integral of the current through SIN junction with gap smearing
	//E = x, V=v, T=tau, [x]=[v]=[tau]
   
   double p, q;
   double a, a1, a2;
   double b, b1;
   double x, dx;
   int l, inf;
   double i; //current
   FILE *f8; //density of states
   
   i=0;
   *dI=0; *dP=0; *dPdI=0;
   double tmpdI=0, tmpdP=0, tmpdPdI=0;
   double tmp1;

   //v=abs(v);

   dx=1e-3;
   inf=(int)(20/dx);

   //f8=fopen("DOS.txt","w");
  
   x=dx;
   #pragma omp parallel for default(shared) private(x, a, b, tmp1, p, q, b1, a2) reduction(+:tmpdI, tmpdP, tmpdPdI) num_threads(THREADS)
	   for (l=0; l<=inf; l++)
	   {  
		   x=dx*(l+1);
		   a=x*x-gamma*gamma-1;	   b=-2*x*gamma;
			tmp1=sqrt(a*a+b*b);

		   p=sqrt((a+tmp1)/2); q=sqrt((tmp1-a)/2);

		   b1=(x*p+gamma*q)/(p*p + q*q); //fprintf(f8,"%f %f\n",x,b1);

		   a1=exp(x/tau); a2=exp((x-v)/tauE); 
		   a=b1*(a1+a2)/(1+a1)/(1+a2);
		   tmpdI+=a; 
		   tmpdP+=x*x*a;
		   tmpdPdI+=-x*a; //fprintf(f8,"%f %f\n",x,a);
	   }

   //fprintf(f8,"\n");
   /*x=-dx;
   for (l=0; l<=inf; l++)
   {  
	   a=x*x-gamma*gamma-1;	   b=-2*x*gamma;

	   p=sqrt((a+sqrt(a*a+b*b))/2); q=-sqrt((sqrt(a*a+b*b)-a)/2);

	   b1=-(x*p+gamma*q)/(p*p + q*q); //fprintf(f8,"%f %f\n",x,b1);

	   a1=exp(x/tau); a2=exp((x-v)/tauE); 
	   a=b1*(a1+a2)/(1+a1)/(1+a2);

	   *dI+=a; //fprintf(f8,"%f %f\n",x,i);
       *dP+=x*x*a;
       *dPdI+=-x*a; fprintf(f8,"%f %f\n",x,a);
	   x-=dx;
   }*/
   //fclose(f8);

   *dI=tmpdI*2*dx;
   *dP=tmpdP*2*dx;
   *dPdI=tmpdPdI*2*dx;
   
   return i;

}


int main()
{
  FILE *f3;  //IV 
  int k,q,j,l,n;
  int Nt, NV;

  float Vol; //volume of absorber in um^3

  float Z; // heat exchange in normal metal in nW/(K^5*micron^3)
 
  float Delta;

  float Rn; //normal resistance of 1 SIN, kOhm 

  float Rsg, Rsg1; //subgap resistance of 1 SIN, kOhm 

  float Rleak; //leakage resistance of SIN

  float Ra; //resistance of absorber, kOhm

  float Pbg; //background power, pW

  float I[10002], I1[10002];
  float V[10002];
  float Ib; //bias for electrone cooling

  double I0; //units of current

  float dV, Vfin; //range of voltage

  float E, dE, E1; //energy, eV

  double int1, int2; //integrals

  float Te, Tp, dT, Ts, Tc; //electron and phonon temperature

  double tau, Vg, tauE, tauC, T1, T2; //dimentionless variables
  double tauold;

  double Pe_p, Pabs, Pleak, Pheat, Pcool, P, Pcool1, Ps;

  int M, MP; //number of bolometers in series, MP - in parallel

  float Rn1, Rleak1, V1, V2, Vcur; //for thermometer junctions

  float gamma, beta; //gap smearing, returned power

  float G, dPdT, dIdT, dIdV, dPdV, Sv, dPT, mm; //for noise 
  double dP, dI, dPdI;
  float NEP, NEPs, NEPep, NEPa, NEPph;
  float vn, in; //amplifyer noise

  float Tp0, x, DeltaT;
  float eps; //coefficient between Ps and Ts

 double ZS, VolS; //sigma and volume of supercond-r

 float NoiA; //amplifier noise

 double dPbgPB=1.288492, dBeta=0.370992, dOddParam=0.040992;  

  clock_t start, finish;
  start=clock();

  f3=fopen("Te.txt","w"); 
  fclose(f3);

      
    // --------- known/guessed physical parameters
   
     M=6;  //in series
	 MP=1;  //in parallel

     Vol = 0.02; VolS = 2.5; //um^3  //0.02 - 1bolo in Ol array, 0.0212 - Sigma (T8)

	 Z = 1.25; ZS = 0.3; //Sigma   nW/K^5/um^3

	 Tc=1.47; 

	 Delta = 1.764*Tc; //in K , Vg[eV]=Tc * 1.764 * 86.25e-6

	 gamma=5e-5; //smearing, check integration step below 1e-4!!!

	 beta=dBeta;   //returning power //<1
	 
     Tp = 0.300; 	//0.31 temperature, K

	 Te = 0.300; //is to be found

	 Ts = 0.300; //electrone T in SC 

	 Rn = 8.05e3/M*MP; //68e3; /*0.5*4.8e3;*/ per (2SINs + abs) 

	 Ra = 112.5;  //Ohm  225/2 - L27, 220/2 - L43

	 Ib = 1.4; //nA, for electron cooling

	 Pbg = dPbgPB/M/MP; //pW, per 1 bolometer //1.4

	// Rn1=57e3/M; //0.5*3.57e3;
	// Rleak1 = Rn1/gamma;

	 //---------- normalized constants

	 Rn=(Rn-Ra)/2;
	 //Rleak=Rn/gamma; //Ohm, is not needed if Gamma<>0

	 I0 = Delta/Rn*KBOL*1e9;

	 Vg = Delta*KBOL;  printf("%f\n",Vg);

	 tau = Ts/Delta;
	 tauE = Te/Delta;
	 tauC =Tc/Delta;

	 //---------- calculation parameters


	 dV=5e-6; //voltage step, V
	 Vfin=0.8*Vg;
	 NV=1*Vfin/dV;
	 for (j=0; j<=NV; j++) V[j]=(j*dV); // = [V]

	 Tp0=Tp; DeltaT=1;

     //----------------------------------------------------------
	 //---------- find Te for the whole IV-curve ----------------

	 f3=fopen("Te.txt","a");
	 FILE* check=fopen("conv.txt", "w");

	 //j=5;
	 for (j=1; j<=NV-1; j++)   //next voltage
	 {
	   
		 for (n=1; n<=5; n++)   //next interation
	     {
		 
	        //-- solve HBE for Absorber ------------------
	   
	        T1=0; T2=3/1.764;
			tauE=(T1+T2)/2;
			DeltaT = sqrt(1-pow(tau/tauC, double(3.2)));
			tauold=0;			

				while(abs(tauold-tauE)/tauE>0.001) //next iteration
				{
				    tauold=tauE;
					//fprintf(check, "%d %lf\n", l, tauE);
					I[j] = 1e-12*currentIntCom(DeltaT,V[j]/Vg,tau,tauE,gamma)*I0; //in nA

					Pe_p = Z * Vol * (pow(Tp,5)-pow(tauE*Delta,5)) * 1e3; //pW

					Pabs = I[j]*I[j] * Ra * 1e-6;   //pW

					//Rsg = 2 * dV/abs(I[j+1]-I[j-1])*1e9; //in Ohm //2 from derivative

					Pcool =  1*PowerCoolIntCom(DeltaT, V[j]/Vg, tau, tauE,gamma, &Ps)* Vg*Vg/Rn * 1e12;  //pW, 0.1 - experimental parameter

					Ps=Ps* Vg*Vg/Rn * 1e12; //returning power from S to N

					Pheat = Pe_p + Pabs + Pbg + 2*beta*Ps;

					P=Pheat - 2*Pcool;
					if (P<0) {T2=tauE;}	else {T1=tauE;}
					tauE=(T1+T2)/2;
				} 

				//-- solve HBE for Superconductor -----------

				T1=0; T2=2/1.764;
			    tau=(T1+T2)/2;

				while(abs(tauold-tau)/tau>0.001) //next iteration
				{
				   tauold=tau;
				   //fprintf(check, "%d %lf\n", l, tau);
				   DeltaT = sqrt(1.0-pow(float(tau*1.764),float(3.2)));

					Pe_p = 0.98 * ZS * VolS * (pow(Tp,5)-pow(tau*Delta,5))*exp(-DeltaT/tau) * 1e3; //pW

					Pcool =  1*PowerCoolIntCom(DeltaT, V[j]/Vg, tau, tauE,gamma, &Ps)* Vg*Vg/Rn * 1e12;  //pW, 0.1 - experimental parameter

					Ps=Ps* Vg*Vg/Rn * 1e12; //returning power from S to N

					P =dOddParam*(1-beta)*Ps + Pe_p;
					if (P<0) {T2=tau;}
					else {T1=tau;}
					tau=(T1+T2)/2;
					//Te = tauE*Delta; Ts = tau*Delta;
				} 
			}

	   //Te = tauE*Delta;
	   //Ts = tau*Delta;
	   DeltaT = sqrt(1-pow(tau/tauC, double(3.2)));
	   //DeltaT=1;

	   I[j]=1e-12*currentIntCom(DeltaT, V[j]/Vg,tau,tauE,gamma)*I0;
	   fprintf(f3,"%f %e\n", M*(2*V[j]+I[j]*Ra*1e-9), I[j]*MP*1e-09); //!!!extra 1e-09 in i
	   printf("%f %e\n", M*(2*V[j]+I[j]*Ra*1e-9), I[j]*MP*1e-09);	//!!!extra 1e-09 in i
	   //fprintf(f4,"%f %f\n", M*(2*V[j]+I[j]*Ra*1e-9), 2*Rsg);  

   
	 }
    
	 fclose(f3);
	  
	   finish=clock();
       double sTime=(double)(finish-start)/CLOCKS_PER_SEC;
       printf("%f\n",sTime);

	getch();
        

}

#undef PI
#undef ME  
#undef E  
#undef HPLANCK
#undef HBAR
#undef KBOL