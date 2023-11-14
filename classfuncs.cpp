#include "header.h"
double CFoo::operator()(double dParam)
{
	double dMinimize;
	par[iParNum]=dParam;

	countnum=CEB_2eq_parallel_lite();
	Resample(countnum, countexp, Irex, Vrex, Iexp, Vexp, Inum, Vnum);

	dMinimize=ChiSqDer(countnum, Vnum, Inum, Irex);
	return dMinimize;
} 
CFoo::CFoo(int parnum)
{/*Used to set the base parameters for fitting*/
	iParNum=parnum;
	struct parametername
	{
		char parname[20];
	};
	FILE *parfile;
	float param[21]; 
	int Fit[21];
	struct parametername parameters[21];
	int i=0;
	parfile = fopen("startparams.txt", "r");
	if (!parfile)
	{
		printf("No parameters file!");
		getchar();
		exit(0);
	}
	while (fscanf (parfile, "%s %f %i", parameters[i].parname, &param[i], &Fit[i]) != EOF) 
	{
		printf("%s = %f, to fit = %i\n", parameters[i].parname, param[i], Fit[i]); 
		i++;
	}
	fclose(parfile);

	par[0]=param[0]; Pbg=par[0]; ToFit[0]=Fit[0];		//power for all structure

	par[1]=param[1]; beta=par[1]; ToFit[1]=Fit[1];		//returned power

	par[2]=param[2]; TephPOW=par[2]; ToFit[2]=Fit[2];	//exponent for Te-ph

	par[3]=param[3]; gamma=par[3]; ToFit[3]=Fit[3];		//gap smearing

	par[4]=param[4]; Vol=par[4]; ToFit[4]=Fit[4];		//volume of absorber

	par[5]=param[5]; VolS=par[5]; ToFit[5]=Fit[5];		//volume of superconductor

	par[6]=param[6]; Z=par[6]; ToFit[6]=Fit[6];		//heat exchange in normal metal

	par[7]=param[7]; ZS=par[7]; ToFit[7]=Fit[7];		//sigma of superconductor

	par[8]=param[8]; Tc=par[8]; ToFit[8]=Fit[8];		//critical temperature of superconductor

	par[9]=param[9]; Rn=par[9]; ToFit[9]=Fit[9];		//total normal resistance

	par[10]=param[10]; Rleak=par[10]; ToFit[10]=Fit[10];	//leakage resistance

	par[11]=param[11]; Wt=par[11]; ToFit[11]=Fit[11];	//transparency of barter

	par[12]=param[12]; tm=par[12]; ToFit[12]=Fit[12];	//depairing energy

	par[13]=param[13]; ii=par[13]; ToFit[13]=Fit[13];	//coefficient for Andreev current

	par[14]=param[14]; Ra=par[14]; ToFit[14]=Fit[14];	//resistance of absorber
	
	par[15]=param[15]; M=par[15]; ToFit[15]=Fit[15];	//number of bolometers in series

	par[16]=param[16]; MP=par[16]; ToFit[16]=Fit[16];	//number of bolometers in parallel

	par[17]=param[17]; Tp=par[17]; ToFit[17]=Fit[17];	//phonon temperature

	par[18]=param[18]; dVFinVg=par[18]; ToFit[18]=Fit[18];
	
	par[19]=param[19]; dVStartVg=par[19]; ToFit[19]=Fit[19];
	
	par[20]=param[20]; dV=par[20]; ToFit[20]=Fit[20];

	char fname[]="SINS1_53_240_303_mK.txt";			//default filename here. Format: two columns (voltage (V), current(A)), decimal delimiter: .
	
	//printf("Enter the name of file: \n");			//example: "ol g7-25nn\200mk_2g.txt"
	//scanf("%[^\n]", fname);				//your filename here	Iexp=NULL;
	
	Vexp=NULL;
	Inum=NULL;
	Vnum=NULL;
	countexp=GetExp(fname, 1, 1, Iexp, Vexp, 0);
	countnum=CEB_2eq_parallel_lite();
	Irex=new double[countnum];
	Vrex=new double[countnum];
	Resample(countnum, countexp, Irex, Vrex, Iexp, Vexp, Inum, Vnum);
	
	FILE* conv=fopen("converg.txt", "w");
	fprintf(conv, "%lf %e 0\n", par[iParNum], ChiSq(countnum, Vnum, Inum, Irex));
	fclose(conv);	
}
	
CFoo::~CFoo()
{
	if (Iexp!=NULL) delete[] Iexp;
	if (Vexp!=NULL) delete[] Vexp;
	if (Inum!=NULL) delete[] Inum;
	if (Vnum!=NULL) delete[] Vnum;
	if (Irex!=NULL) delete[] Irex;
	if (Vrex!=NULL) delete[] Vrex;
}

void CFoo::SeqFit(int iRunCount)
{/*Fits using Golden method*/
	Golden method(1e-3);

	for (int i=0; i<iRunCount; i++)
	{
		int ParSeq[iNumParams];
		int Tmp[iNumParams];
		int iRandom=0;

		srand(time(NULL));
		for (int j=0; j<iNumParams; j++) Tmp[j]=j;
		for (int j=0; j<iNumParams; j++)
		{
			iRandom=rand()%(iNumParams-j);
			ParSeq[j]=Tmp[iRandom];
			for (int k=0; k<iNumParams-iRandom-1; k++) 
			{
				Tmp[iRandom+k]=Tmp[iRandom+k+1];
			}
		}
		//for (int j=0; j<iNumParams; j++) ParSeq[j]=j;

		for (int j=0; j<iNumParams; j++) if (ToFit[ParSeq[j]])
		{
			iParNum=ParSeq[j];
			method.ax=0.5*par[iParNum];
			method.bx=1*par[iParNum];
			method.cx=2*par[iParNum];
			//printf("Bracketed at: %lf %lf %lf", brent.ax, brent.bx, brent.cx);
			//scanf("%lf", &tmp);
			
			par[iParNum]=method.minimize(*this);
		}
		FILE* params=fopen("fitparameters_new.txt", "a");
		for (int i=0; i<iNumParams; i++) fprintf(params, "%lf\t", par[i]);
		fprintf(params, "%e\n", method.fmin);
		fclose(params);
	}
}
COff::COff(double *Iexp, double *Vexp, long countexp)
{
	count=countexp;
	Iofx=new double[countexp];
	Vofx=new double[countexp];
	Iref=new double[countexp];
	Vref=new double[countexp];
	for (long i=0; i<countexp; i++)
	{
		Vref[i]=Vexp[i];
		Iref[i]=Iexp[i];
	}
}
COff::~COff()
{
	if (Iofx!=NULL) delete[] Iofx;
	if (Vofx!=NULL) delete[] Vofx;
	if (Iref!=NULL) delete[] Iref;
	if (Vref!=NULL) delete[] Vref;
}
double COff::operator ()(double dOffset)
{
	long lLowI=0;
	double tmp=fabs(Iref[0]);
	for (long i=0; i<count; i++)
	{
		Vofx[i]=Vref[i]-dOffset;
		Iofx[i]=Iref[i];
		if (fabs(Iref[i])<tmp)
		{
			lLowI=i;
			tmp=fabs(Iref[i]);
		}
	}
	if (Iofx[lLowI]<0) lLowI++;
	long lLengthPos=count-lLowI;
	long lLengthNeg=count-lLengthPos;
	double *Vlow, *Ilow, *Vhigh, *Ihigh;
	long lLengthLow=0, lLengthHigh=0;
	if (lLengthPos>lLengthNeg)
	{
		lLengthHigh=lLengthPos;
		lLengthLow=lLengthNeg;
		Vhigh=new double[lLengthHigh];
		Ihigh=new double[lLengthHigh];
		Vlow=new double[lLengthLow];
		Ilow=new double[lLengthLow];
		for (long i=0; i<lLengthHigh; i++) 
		{
			Vhigh[i]=fabs(Vofx[lLowI+i]);
			Ihigh[i]=fabs(Iofx[lLowI+i]);
		}
		for (long i=0; i<lLengthLow; i++) 
		{
			Vlow[i]=fabs(Vofx[lLowI-1-i]);
			Ilow[i]=fabs(Iofx[lLowI-1-i]);
		}		
	}
	else
	{
		lLengthHigh=lLengthNeg;
		lLengthLow=lLengthPos;
		Vhigh=new double[lLengthHigh];
		Ihigh=new double[lLengthHigh];
		Vlow=new double[lLengthLow];
		Ilow=new double[lLengthLow];
		for (long i=0; i<lLengthHigh; i++) 
		{
			Vhigh[i]=fabs(Vofx[lLowI-1-i]);
			Ihigh[i]=fabs(Iofx[lLowI-1-i]);
		}
		for (long i=0; i<lLengthLow; i++) 
		{
			Vlow[i]=fabs(Vofx[lLowI+i]);
			Ilow[i]=fabs(Iofx[lLowI+i]);
		}	
	}
	double *Irex=new double[lLengthLow], *Vrex=new double[lLengthLow];
	FILE* test1=fopen("test1.txt", "w+");
	FILE* test2=fopen("test2.txt", "w+");
	FILE* test3=fopen("test3.txt", "w+");
	for (int i=0; i<lLengthLow; i++) fprintf(test1, "%e %e \n", Vlow[i], Ilow[i]);
	for (int i=0; i<lLengthHigh; i++) fprintf(test2, "%e %e \n", Vhigh[i], Ihigh[i]);
	Resample(lLengthLow, lLengthHigh, Vrex, Irex, Vhigh, Ihigh, Vlow, Ilow);
	for (int i=0; i<lLengthLow; i++) fprintf(test3, "%e %e \n", Vrex[i], Irex[i]);
	fclose(test1);
	fclose(test2);
	fclose(test3);
	double dResult=ChiSqHi(lLengthLow, Ilow, Vlow, Vrex);
	if (Ilow!=NULL) delete[] Ilow;
	if (Vlow!=NULL) delete[] Vlow;
	if (Ihigh!=NULL) delete[] Ihigh;
	if (Vhigh!=NULL) delete[] Vhigh;
	if (Irex!=NULL) delete[] Irex;
	if (Vrex!=NULL) delete[] Vrex;
	return dResult;
}
