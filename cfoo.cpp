#include "cfoo.h"
#include "CEB_2eq_parallel_lite.h"
#include "constants.h"
#include "mins.h"
#include "utils.h"

double CFoo::operator()(double dParam)
{
    double dMinimize;
    par[iParNum] = dParam;

    countnum = CEB_2eq_parallel_lite();
    Resample(countnum, countexp, Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    dMinimize = ChiSqDer(countnum, Vnum, Inum, Irex);
    return dMinimize;
}

CFoo::CFoo(int parnum)
{ /*Used to set the base parameters for fitting*/
    iParNum = parnum;
    struct parametername
    {
        char parname[20];
    };
    FILE *parfile;
    float param[21];
    int Fit[21];
    struct parametername parameters[21];
    int i = 0;
    parfile = fopen("startparams.txt", "r");
    if (!parfile) {
        printf("No parameters file!");
        getchar();
        exit(0);
    }
    while (fscanf(parfile, "%s %f %i", parameters[i].parname, &param[i], &Fit[i]) != EOF) {
        printf("%s = %f, to fit = %i\n", parameters[i].parname, param[i], Fit[i]);
        i++;
    }
    fclose(parfile);

    par[0] = param[0];
    Pbg = par[0];
    ToFit[0] = Fit[0]; //power for all structure

    par[1] = param[1];
    beta = par[1];
    ToFit[1] = Fit[1]; //returned power

    par[2] = param[2];
    TephPOW = par[2];
    ToFit[2] = Fit[2]; //exponent for Te-ph

    par[3] = param[3];
    gamma = par[3];
    ToFit[3] = Fit[3]; //gap smearing

    par[4] = param[4];
    Vol = par[4];
    ToFit[4] = Fit[4]; //volume of absorber

    par[5] = param[5];
    VolS = par[5];
    ToFit[5] = Fit[5]; //volume of superconductor

    par[6] = param[6];
    Z = par[6];
    ToFit[6] = Fit[6]; //heat exchange in normal metal

    par[7] = param[7];
    ZS = par[7];
    ToFit[7] = Fit[7]; //sigma of superconductor

    par[8] = param[8];
    Tc = par[8];
    ToFit[8] = Fit[8]; //critical temperature of superconductor

    par[9] = param[9];
    Rn = par[9];
    ToFit[9] = Fit[9]; //total normal resistance

    par[10] = param[10];
    Rleak = par[10];
    ToFit[10] = Fit[10]; //leakage resistance

    par[11] = param[11];
    Wt = par[11];
    ToFit[11] = Fit[11]; //transparency of barter

    par[12] = param[12];
    tm = par[12];
    ToFit[12] = Fit[12]; //depairing energy

    par[13] = param[13];
    ii = par[13];
    ToFit[13] = Fit[13]; //coefficient for Andreev current

    par[14] = param[14];
    Ra = par[14];
    ToFit[14] = Fit[14]; //resistance of absorber

    par[15] = param[15];
    M = par[15];
    ToFit[15] = Fit[15]; //number of bolometers in series

    par[16] = param[16];
    MP = par[16];
    ToFit[16] = Fit[16]; //number of bolometers in parallel

    par[17] = param[17];
    Tp = par[17];
    ToFit[17] = Fit[17]; //phonon temperature

    par[18] = param[18];
    dVFinVg = par[18];
    ToFit[18] = Fit[18];

    par[19] = param[19];
    dVStartVg = par[19];
    ToFit[19] = Fit[19];

    par[20] = param[20];
    dV = par[20];
    ToFit[20] = Fit[20];

    char fname[]
        = "SINS1_53_240_303_mK.txt"; //default filename here. Format: two columns (voltage (V), current(A)), decimal delimiter: .

    //printf("Enter the name of file: \n");			//example: "ol g7-25nn\200mk_2g.txt"
    //scanf("%[^\n]", fname);				//your filename here	Iexp=NULL;

    Vexp = NULL;
    Inum = NULL;
    Vnum = NULL;
    countexp = GetExp(fname, 1, 1, Iexp, Vexp, 0);
    countnum = CEB_2eq_parallel_lite();
    Irex = new double[countnum];
    Vrex = new double[countnum];
    Resample(countnum, countexp, Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    FILE *conv = fopen("converg.txt", "w");
    fprintf(conv, "%lf %e 0\n", par[iParNum], ChiSq(countnum, Vnum, Inum, Irex));
    fclose(conv);
}

CFoo::~CFoo()
{
    if (Iexp != NULL)
        delete[] Iexp;
    if (Vexp != NULL)
        delete[] Vexp;
    if (Inum != NULL)
        delete[] Inum;
    if (Vnum != NULL)
        delete[] Vnum;
    if (Irex != NULL)
        delete[] Irex;
    if (Vrex != NULL)
        delete[] Vrex;
}

void CFoo::SeqFit(int iRunCount)
{ /*Fits using Golden method*/
    Golden method(1e-3);

    for (int i = 0; i < iRunCount; i++) {
        int ParSeq[iNumParams];
        int Tmp[iNumParams];
        int iRandom = 0;

        srand(time(NULL));
        for (int j = 0; j < iNumParams; j++)
            Tmp[j] = j;
        for (int j = 0; j < iNumParams; j++) {
            iRandom = rand() % (iNumParams - j);
            ParSeq[j] = Tmp[iRandom];
            for (int k = 0; k < iNumParams - iRandom - 1; k++) {
                Tmp[iRandom + k] = Tmp[iRandom + k + 1];
            }
        }
        //for (int j=0; j<iNumParams; j++) ParSeq[j]=j;

        for (int j = 0; j < iNumParams; j++)
            if (ToFit[ParSeq[j]]) {
                iParNum = ParSeq[j];
                method.ax = 0.5 * par[iParNum];
                method.bx = 1 * par[iParNum];
                method.cx = 2 * par[iParNum];
                //printf("Bracketed at: %lf %lf %lf", brent.ax, brent.bx, brent.cx);
                //scanf("%lf", &tmp);

                par[iParNum] = method.minimize(*this);
            }
        FILE *params = fopen("fitparameters_new.txt", "a");
        for (int i = 0; i < iNumParams; i++)
            fprintf(params, "%lf\t", par[i]);
        fprintf(params, "%e\n", method.fmin);
        fclose(params);
    }
}

long CFoo::CEB_2eq_parallel_lite(void)
{
    FILE *f2, *f3, *f3old, *f4, *f5; //Noise, Te, Te_old, NEP, G

    int k, q, j, l, n; //counters

    int Nt, NV; //numbers of temperature and voltage steps

    double I[10002], I1[10002];

    float I_A[10002], I_As[10002];

    double V[10002];

    double Ib; //bias for electron cooling

    double I0; //units of current

    double Vstr, Vfin; //range of voltage

    float Z; //heat exchange in normal metal in nW / (K^5 * micron^3)

    double Delta, gamma; //energy gap and gap smearing

    double E, dE, E1; //energy, eV

    double int1, int2; //integrals

    double Pbg; //background power, pW

    double Te, dT, Ts; //electron and phonon temperature

    double tau, Vg, tauE, tauC, T1, T2; //dimentionless variables

    double tauold;

    double Pe_p, Pabs, Pleak, Pheat, Pcool, P, Pcool1, Ps, Ps1, Pand; //for power

    double Rsg, Rsg1; //subgap resistance of 1 SIN, kOhm

    double Rleak; //leakage resistance of SIN

    double Rsin; //normal resistance of single SIN junction

    double Rn1, Rleak1, V1, V2, Vcur; //for thermometer junctions

    double G, dPdT, dIdT, dIdV, dPdV, Sv, dPT, mm; //for noise

    double dP, dI, dPdI; //for power integrals

    double NEP, NEPs, NEPep, NEPa, NEPph; //noise equivalent power

    double vn, in; //amplifier voltage and current noise

    double Tp0, x, DeltaT;

    double eps; //coefficient between Ps and Ts

    double NoiA; //amplifier noise

    float G_NIS, G_e; //thermal conductivity

    float Wt; //transparency of barrier

    float tm; //depairing energy

    double I_A0, I_As0;

    float ii; //coefficient for Andreev current

    clock_t start, finish;

    start = clock();

    char c = 0;

    // --------- known/guessed physical parameters

    Pbg = par[0]; //incoming power for all structure, pW

    beta = par[1]; //returning power ratio, <1

    TephPOW = par[2]; //exponent for Te-ph

    // gamma = par[3]; //gap smearing, not used, should be 0

    Vol = par[4]; //volume of absorber, um^3

    VolS = par[5]; //volume of superconductor, um^3

    Z = par[6]; //heat exchange in normal metal, nW/(K^5*um^3)

    ZS = par[7]; //sigma

    Tc = par[8]; //critical temperature, K

    Rn = par[9] / M * MP; //normal resistance for 1 bolometer, Ohm

    Rleak = par[10] / M * MP; //leakage resistance for 1 bolometer, Ohm

    Wt = par[11]; //transparency of barrier

    tm = par[12]; //depairing energy

    ii = par[13]; //coefficient for Andreev current

    Ra = par[14]; //normal resistance for 1 absorber, Ohm

    M = par[15]; //number of bolometers in series

    MP = par[16]; //number of bolometers in parallel

    Tp = par[17]; //phonon temperature, K

    dVFinVg = par[18]; //voltage range end

    dVStartVg = par[19]; //voltage range start

    dV = par[20]; //voltage step, V

    Te = Tp; //electron temperature, is to be found, K

    Ts = Tp; //electron temperature in superconductor, K

    //Ib = 1.4;								//bias current for electron cooling, nA

    dPbg = Pbg / M / MP; //incoming power for 1 bolometer, pW

    Delta = 1.764 * Tc; //in K, Vg[eV] = Tc * 1.764 * 86.25e-6

    f3old = fopen("Te_old.txt", "w");
    f3 = fopen("Te.txt", "r");
    //f3 = fopen(cPbgTe, "r");

    if (f3 != NULL) {
        while (fscanf(f3, "%c", &c) != EOF)
            fprintf(f3old, "%c", c);
        fclose(f3);
    }

    fclose(f3old);

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
    //Rleak=Rn/gamma;							//Ohm, is not needed if Gamma<>0

    I0 = Delta / Rn * KBOL * 1e9;

    Vg = Delta / 11604.505;
    /*Vg = Delta * KBOL;*/

    printf("%lf\n", Vg);

    tau = Ts / Delta;

    tauE = Te / Delta;

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
        V[j] = Vstr + ((j - 0) * dV); // = [V]

    f2 = fopen("Noise.dat", "a");
    f3 = fopen("Te.txt", "a");
    f4 = fopen("NEP.dat", "a");
    f5 = fopen("G.txt", "a");

    fprintf(f2, "Voltage\tNOISEep\tNOISEs\tNOISEa\tNOISE\tNOISEph\tNOISE^2-NOISEph^2\n");
    fprintf(f3, "Voltage\tCurrent\tIqp\tIand\tV/Rleak\tTe\tTs\tDeltaT\tPand\tPleak\tPabs\tPcool\n");
    fprintf(f4, "Voltage\tCurrent\tNEPep\tNEPs\tNEPa\tNEP\tNEPph\tSv\tNEP^2-NEPph^2\n");
    fprintf(f5, "Voltage\tGe\tGnis\n");

    for (j = 1; j <= NV - 1; j++) //next voltage
    {
        for (n = 1; n <= 5; n++) //next interation
        {
            T1 = 0;

            T2 = 3 / 1.764;

            tauE = (T1 + T2) / 2;

            DeltaT = std::sqrt(1 - pow(Ts / Tc, float(3.2)));

            for (l = 1; l <= 15; l++) //next iteration
            {
                tauE = (T1 + T2) / 2;

                I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9; //in nA

                I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0; //in nA

                Pe_p = Z * Vol
                       * (std::pow(Tp, TephPOW /*7,6,5*/)
                          - std::pow(tauE * Delta, TephPOW /*7,6,5*/))
                       * 1e3; //pW

                Pabs = I[j] * I[j] * Ra * 1e-6; //pW

                Pleak = 2.0 * V[j] * V[j] / Rleak * 1e12; // pW, 2 because of 2 SIN

                Pand = (I_A[j] * 1e-3) /*uA*/ * (I_A[j] * 1e-3) /*uA*/ * Ra /*Ohm*/
                       + 2.0 * (I_A[j] * 1e3) /*pA*/ * V[j] /*V*/;          //pW, absorber + Andreev

                Pcool = 1.0 * PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE, &Ps) * Vg * Vg / Rn
                        * 1e12; //pW, 0.1 - experimental parameter

                Ps = Ps * Vg * Vg / Rn * 1e12; //returning power from S to N

                Pheat = Pe_p + 1.0 * Pabs + 1.0 * Pand + dPbg + 2.0 * beta * Ps
                        + 1.0 * Pleak; // + Pand;

                P = Pheat - 2.0 * Pcool;

                if (P < 0) {
                    T2 = tauE;
                }

                else {
                    T1 = tauE;
                }
            }
        }
        //DeltaT = 1; // !!!

        Te = tauE * Delta;

        Ts = tau * Delta;

        DeltaT = std::sqrt(1 - pow(Ts / Tc, float(3.2)));

        I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9;

        I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0;

        fprintf(f3,
                "%f\t%g\t%g\t%g\t%f\t%f\t%g\t%g\t%g\t%g\t%g\t%g\n",
                M * (2 * V[j] + (I[j] * 1e-9 + 1.0e-9 * I_A[j]) * Ra),
                MP * 1e-9 * (I[j] + I_A[j]),
                I[j] * MP * 1e-9,
                I_A[j] * MP * 1e-9,
                (V[j] / Rleak * 1e9) * MP,
                Te,
                Ts,
                DeltaT,
                Pand,
                Pleak,
                Pabs,
                Pcool);

        printf("Voltage: %g Current: %g\n",
               M * (2 * V[j] + (I[j] * 1e-9 + 1.0e-9 * I_A[j]) * Ra),
               MP * 1e-9 * (I[j] + I_A[j]) /*, I[j] * MP * 1e-9, I_A[j] * MP * 1e-9, Pand, Pleak*/);

        //fprintf(f4,"%f %f\n", M * (2 * V[j] + I[j] * Ra * 1e-9), 2 * Rsg);

        Inum[j - 1] = MP * 1e-9 * (I[j] + I_A[j]);

        Vnum[j - 1] = M * (2 * V[j] + ((I[j] + I_A[j]) * Ra * 1e-9));

        //----- NEP ----------------------------------------------

        dT = 0.005;

        vn = (3.2e-9) * std::sqrt(2.0);      //V/sqrt(Hz), for 2 amp * sqrt(2) AD745 amplifier
        in = (6.9e-15) * 1 / std::sqrt(2.0); //A/sqrt(Hz), for 2 amp * 1 / sqrt(2)

        //		vn = (8.0e-9) * sqrt(2.0);					//V/sqrt(Hz), for 2 amp * sqrt(2) OPA111 amplifier
        //		in = (0.8e-15) * 1 / sqrt(2.0);					//A/sqrt(Hz), for 2 amp * 1 / sqrt(2)

        //		vn = (0.9e-9) * sqrt(2.0);					//V/sqrt(Hz), for 2 amp * sqrt(2) AD797 amplifier
        //		in = (2.0e-12) * 1 / sqrt(2.0);					//A/sqrt(Hz), for 2 amp * 1 / sqrt(2)

        //		vn = (1.1e-9) * sqrt(2.0);					//V/sqrt(Hz), for 2 amp * sqrt(2) IFN146 amplifier
        //		in = (0.3e-15) * 1 / sqrt(2.0);					//A/sqrt(Hz), for 2 amp * 1 / sqrt(2)

        //		vn = (5.1e-9) * sqrt(2.0);					//V/sqrt(Hz), for 2 amp * sqrt(2) OPA1641 amplifier
        //		in = (0.8e-15) * 1 / sqrt(2.0);					//A/sqrt(Hz), for 2 amp * 1 / sqrt(2)

        dPT = PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta, &Ps)
              - PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta, &Ps);

        dPdT = Vg * Vg / Rn * 1e12 * (dPT) / (2 * dT); //pW/K

        dIdT = I0
               * (currentInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta)
                  - currentInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta))
               / (2 * dT); //nA/K

        dIdV = (I0
                * (currentInt(DeltaT, V[j + 1] / Vg, tau, tauE)
                   + ii * AndCurrent(DeltaT, V[j + 1] / Vg, tauE, Wt, tm)
                   - currentInt(DeltaT, V[j - 1] / Vg, tau, tauE)
                   - ii * AndCurrent(DeltaT, V[j - 1] / Vg, tauE, Wt, tm)))
               / (2 * dV); //nA/V

        dPdV = Vg * Vg / Rn * 1e12
               * (PowerCoolInt(DeltaT, V[j + 1] / Vg, tau, tauE, &Ps)
                  - PowerCoolInt(DeltaT, V[j - 1] / Vg, tau, tauE, &Ps))
               / (2 * dV); //pW/V

        G_NIS = dPdT;
        G_e = 5 * Z * Vol * std::pow(Te, 4) * 1e3; //pW/K

        G = G_e + 2 * (G_NIS - dIdT / dIdV * dPdV); //pW/K, '2' for 2 SINs

        Sv = -2 * dIdT / dIdV / G; //V/pW, for 1 bolo

        Sv = Sv / MP;

        NEPep = 10 * EL * KBOL * Z * Vol * (std::pow(Tp, TephPOW) + std::pow(Te, TephPOW)) * 1e3
                * 1e12; //^2 pW^2/Hz

        NoiA = (vn * vn + std::pow(in * (2 * 1e9 / dIdV + Ra) * M / MP, 2)); //(V/sqrt(Hz))^2

        NEPa = (NoiA) / Sv / Sv; //pW^2/Hz

        //----- NEP SIN approximation ----------------------------

        dI = 2 * EL * std::abs(I[j]) / std::pow(dIdV * Sv, 2) * 1e9; //pW^2/Hz

        dPdI = 2 * 2 * EL * Pcool / (dIdV * Sv)
               * 1e9; //pW^2/Hz, second '2' is from comparison with integral

        mm = std::log(std::sqrt(2 * PI * KBOL * Te * Vg) / (2 * std::abs(I[j]) * Rn * 1e-9));

        dP = (0.5 + mm * mm) * (KBOL * Te) * (KBOL * Te) * std::abs(I[j]) * EL * 1e-9
             * 1e24; //pW^2/Hz

        NEPs = 2 * (dI - 2 * dPdI + dP); //pW^2/Hz, '2' for 2 SINs, all terms positive

        //fprintf(f5, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);

        //----- NEP SIN integral ---------------------------------
        /* 
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / std::pow(dIdV * Sv, 2) * 1e9;			//pW^2/Hz
  
		dPdI = dPdI / (dIdV * Sv) * EL * Vg * Vg / Rn * 1e12 * 1e9;	//pW^2/Hz
  
		dP = dP * Vg * Vg * Vg / Rn * EL * 1e24;			//pW^2/Hz 
        
		NEPs = 2 * (dI - 2 * dPdI + dP);				//pW^2/Hz, '2' for 2 SINs, all terms positive
  
		fprintf(f9, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);
*/
        //--------------------------------------------------------

        NEPph = std::sqrt(M * MP * 2 * Pbg * 0 * 1e9 * HPLANCK * 1e-12
                          + std::pow(M * MP * Pbg * 1e-12, 2) / 1e3 / 1e9); //W/sqrt(Hz), at 0 GHz

        //NEPph = sqrt(M * MP * 2 * Pbg * 350 * 1e9 * HPLANCK * 1e-12 + pow(M * MP * Pbg * 1e-12, 2) / 1.552 / 1e9);	//W/sqrt(Hz), at 350 GHz

        NEPph = NEPph * 1e12; //pW

        NEP = std::sqrt(M * MP * (NEPep + NEPs) + NEPa + NEPph * NEPph); //all squares

        fprintf(f2,
                "%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                M * (2 * V[j] + I[j] * Ra * 1e-9),
                std::sqrt(M * MP * NEPep) * 1e9 * std::abs(Sv),
                std::sqrt(M * MP * NEPs) * 1e9 * std::abs(Sv),
                std::sqrt(NoiA) * 1e9,
                NEP * 1e9 * std::abs(Sv),
                NEPph * 1e9 * std::abs(Sv),
                1e9 * std::abs(Sv) * std::sqrt(NEP * NEP - NEPph * NEPph));

        fprintf(f4,
                "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                M * (2 * V[j] + I[j] * Ra * 1e-9),
                MP * 1e-9 * (I[j]),
                std::sqrt(M * MP * NEPep) * 1e-12,
                std::sqrt(M * MP * NEPs) * 1e-12,
                std::sqrt(NEPa) * 1e-12,
                NEP * 1e-12,
                NEPph * 1e-12,
                std::abs(Sv) * 1e12,
                1e-12 * std::sqrt(NEP * NEP - NEPph * NEPph));

        fprintf(f5, "%g\t%g\t%g\n", M * (2 * V[j] + I[j] * Ra * 1e-9), G_e, G_NIS);

        printf("Sv: %g Te: %g NEPs: %g NEPt: %g\n\n",
               std::abs(Sv) * 1e12,
               Te,
               std::sqrt(M * MP * NEPs) * 1e-12,
               NEP * 1e-12);

        //system("cls");
    }

    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    finish = clock();
    double sTime = (double) (finish - start) / CLOCKS_PER_SEC;
    printf("\nTime spent: %f sec.", sTime);

    return NV - 1;
}
