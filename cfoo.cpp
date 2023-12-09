#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "CEB_2eq_parallel_lite.h"
#include "constants.h"
#include "mins.h"
#include "utils.h"

#include "cfoo.h"

CFoo::CFoo(int parnum)
{
    /*
     * Set the base parameters for fitting
     */
    par.resize(iNumParams);
    ToFit.resize(iNumParams);

    iParNum = parnum;

    try {
        std::ifstream parfile;
        parfile.open("startparams.txt", std::ifstream::in);

        std::string parname;
        for (int i = 0; i < iNumParams; ++i) {
            parfile >> std::skipws >> parname >> par[i] >> ToFit[i];
            std::clog << parname << " = " << par[i] << ", to fit = " << std::boolalpha << ToFit[i]
                      << std::endl;
        }
        parfile.close();
    } catch (std::exception &e) {
        std::cerr << "Error while reading parameters: " << e.what() << std::endl;
        getchar();
        exit(0);
    }

    /*
     * Default filename here.
     * Format: two columns (voltage (V), current(A)).
     * Decimal delimiter: .
     */
    std::string fname = "SINS1_53_240_303_mK.txt";

    //printf("Enter the name of file: \n");			//example: "ol g7-25nn\200mk_2g.txt"
    //scanf("%[^\n]", fname);				//your filename here	Iexp=NULL;

    GetExp(fname, Iexp, Vexp, false);
    size_t countnum = CEB_2eq_parallel_lite();
    Irex.resize(countnum);
    Vrex.resize(countnum);

    Resample(Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    std::ofstream conv;
    conv.open("converg.txt", std::ofstream::out);
    conv << par[iParNum] << ChiSq(Inum, Irex) << 0.0 << std::endl;
    conv.close();
}

double CFoo::operator()(double dParam)
{
    par[iParNum] = dParam;

    CEB_2eq_parallel_lite();
    Resample(Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    return ChiSqDer(Vnum, Inum, Irex);
}

void CFoo::SeqFit(int iRunCount)
{
    /*
     * Fit using Golden method
     */

    Golden method(1e-3);

    std::default_random_engine generator;

    for (int i = 0; i < iRunCount; i++) {
        std::vector<int> ParSeq(iNumParams);
        std::vector<int> Tmp(iNumParams);

        // fill `Tmp` with a sequence of ints from 0 to (iNumParams - 1)
        std::iota(Tmp.begin(), Tmp.end(), 0);

        srand(time(NULL));
        int iRandom;
        for (int j = 0; j < iNumParams; ++j) {
            std::uniform_int_distribution<int> distribution(0, iNumParams - j);
            iRandom = distribution(generator);
            ParSeq[j] = Tmp[iRandom];
            // shift items of `Tmp` from `iRandom` to the end
            for (int k = iRandom + 1; k < iNumParams; ++k) {
                Tmp[k - 1] = Tmp[k];
            }
        }
        //for (int j=0; j<iNumParams; j++) ParSeq[j]=j;

        for (int j = 0; j < iNumParams; ++j)
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
    double Ib; //bias for electron cooling

    double I0; //units of current

    double Vstr; // initial voltage
    double Vfin; // final voltage

    double E, dE, E1; //energy [eV]

    double int1, int2; //integrals

    double tau, Vg, tauE, tauC, T1, T2; //dimentionless variables

    double tauold;

    double Pe_p, Pabs, Pleak, Pheat, Pcool, P, Pcool1, Ps, Ps1, Pand; //for power

    double Rsg, Rsg1; //subgap resistance of 1 SIN, kOhm

    double Rsin; //normal resistance of single SIN junction

    double Rn1, Rleak1, V1, V2, Vcur; //for thermometer junctions

    double G, dPdT, dIdT, dIdV, dPdV, Sv, dPT, mm; //for noise

    double dP, dI, dPdI; //for power integrals

    double NEP, NEPs, NEPep, NEPa, NEPph; //noise equivalent power

    double Tp0, x, DeltaT;

    double eps; //coefficient between Ps and Ts

    double NoiA; //amplifier noise

    float G_NIS, G_e; //thermal conductivity

    double I_A0, I_As0;

    clock_t start = clock();

    // --------- known/guessed physical parameters

    double Pbg = par[0]; // incoming power for all structure [pW]

    double beta = par[1]; // returning power ratio, <1

    double TephPOW = par[2]; // exponent for Te-ph

    // double gamma = par[3]; // gap smearing, not used, should be 0

    double Vol = par[4]; // volume of the absorber [um³]

    // double VolS = par[5]; // volume of the superconductor [um³]

    double Z = par[6]; // heat exchange in normal metal [nW/(K⁵×um³)]

    // double ZS = par[7]; // sigma

    double Tc = par[8]; // critical temperature [K]

    size_t M = par[15]; // number of bolometers in series

    size_t MP = par[16]; // number of bolometers in parallel

    double Rn = par[9] / M * MP; // normal resistance for 1 bolometer [Ohm]

    double Rleak = par[10] / M * MP; // leakage resistance per 1 bolometer [Ohm]

    double Wt = par[11]; // transparency of the barrier

    double tm = par[12]; // depairing energy

    double ii = par[13]; // coefficient for Andreev current

    double Ra = par[14]; // normal resistance of 1 absorber [Ohm]

    double Tp = par[17]; // phonon temperature [K]

    double dVFinVg = par[18]; // voltage range end [V]

    double dVStartVg = par[19]; // voltage range start [V]

    double dV = par[20]; // voltage step [V]

    double Te = Tp; // electron temperature to be found [K]
    double Ts = Tp; // electron temperature in superconductor [K]
    const double dT = 0.005;

    //double Ib = 1.4;								//bias current for electron cooling [nA]

    double dPbg = Pbg / M / MP; // incoming power per 1 bolometer [pW]

    double Delta = 1.764 * Tc; // energy gap [K], Vg[eV] = Tc * 1.764 * 86.25e-6

    // if there is a file named “Te.txt”, backup its content into “Te_old.txt”
    {
        FILE *file_Te = fopen("Te.txt", "r");
        //FILE *file_Te = fopen(cPbgTe, "r");

        if (file_Te != NULL) {
            // backup the data
            FILE *file_Te_old = fopen("Te_old.txt", "w");
            char c;
            while (fscanf(file_Te, "%c", &c) != EOF)
                fprintf(file_Te_old, "%c", c);
            fclose(file_Te);
            fclose(file_Te_old);
        }
    }

    FILE *file_Noise = fopen("Noise.dat", "w");
    FILE *file_Te = fopen("Te.txt", "w");
    FILE *file_NEP = fopen("NEP.dat", "w");
    FILE *file_G = fopen("G.txt", "w");

    fclose(file_Noise);
    fclose(file_Te);
    fclose(file_NEP);
    fclose(file_G);

    //---------- normalized constants

    Rn = (Rn - Ra) / 2.0;
    // Rleak = Rn / gamma;							//Ohm, is not needed if Gamma<>0

    I0 = Delta / Rn * KBOL * 1e9; // nA

    Vg = Delta * KBOL;

    printf("Vg = %lf\n", Vg);

    tau = Ts / Delta;
    tauE = Te / Delta;

    //---------- calculation parameters

    Vfin = dVFinVg * Vg;
    Vstr = dVStartVg * Vg;

    size_t NV = Vfin / dV; // the numbers of voltage steps

    Inum.resize(NV - 1);
    Vnum.resize(NV - 1);

    std::valarray<double> V(NV + 1); // [V]
    std::iota(begin(V), end(V), 0.0);
    V = Vstr + (V * dV);

    file_Noise = fopen("Noise.dat", "a");
    file_Te = fopen("Te.txt", "a");
    file_NEP = fopen("NEP.dat", "a");
    file_G = fopen("G.txt", "a");

    fprintf(file_Noise, "Voltage\tNOISEep\tNOISEs\tNOISEa\tNOISE\tNOISEph\tNOISE^2-NOISEph^2\n");
    fprintf(file_Te,
            "Voltage\tCurrent\tIqp\tIand\tV/Rleak\tTe\tTs\tDeltaT\tPand\tPleak\tPabs\tPcool\n");
    fprintf(file_NEP, "Voltage\tCurrent\tNEPep\tNEPs\tNEPa\tNEP\tNEPph\tSv\tNEP^2-NEPph^2\n");
    fprintf(file_G, "Voltage\tGe\tGnis\n");

    std::valarray<double> I(NV + 1);
    std::valarray<double> I_A(NV + 1);

    for (size_t j = 1; j < NV; ++j) // next voltage; exclude edges
    {
        for (size_t n = 0; n < 5; ++n) // next interation
        {
            T1 = 0;
            T2 = 3 / 1.764;
            tauE = (T1 + T2) / 2;
            DeltaT = std::sqrt(1 - pow(Ts / Tc, float(3.2)));

            for (size_t l = 0; l < 15; l++) // next iteration
            {
                tauE = (T1 + T2) / 2.0;

                I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9; // nA

                I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0; // nA

                Pe_p = Z * Vol
                       * (std::pow(Tp, TephPOW /*7,6,5*/)
                          - std::pow(tauE * Delta, TephPOW /*7,6,5*/))
                       * 1e3; // pW

                Pabs = std::pow(I[j], 2) * Ra * 1e-6; // pW

                Pleak = 2.0 * std::pow(V[j], 2) / Rleak * 1e12; // pW, 2 for 2 SINs

                Pand = std::pow(I_A[j] * 1e-3 /*uA*/, 2) * Ra      /*Ohm*/
                       + 2.0 * (I_A[j] * 1e3) /*pA*/ * V[j] /*V*/; //pW, absorber + Andreev

                std::tie(Pcool, Ps) = PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE);

                Pcool *= std::pow(Vg, 2) / Rn * 1e12; // pW
                Ps *= std::pow(Vg, 2) / Rn * 1e12;    // returning power from S to N

                Pheat = Pe_p + Pabs + Pand + dPbg + 2.0 * beta * Ps + Pleak; // + Pand;

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

        DeltaT = std::sqrt(1.0 - std::pow(Ts / Tc, 3.2));

        I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9;

        I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0;

        fprintf(file_Te,
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
               M * (2 * V[j] + 1.0e-9 * (I[j] + I_A[j]) * Ra),
               MP * 1e-9 * (I[j] + I_A[j]) /*, I[j] * MP * 1e-9, I_A[j] * MP * 1e-9, Pand, Pleak*/);

        //fprintf(file_NEP,"%f %f\n", M * (2 * V[j] + I[j] * Ra * 1e-9), 2 * Rsg);

        Inum[j - 1] = MP * 1e-9 * (I[j] + I_A[j]);
        Vnum[j - 1] = M * (2.0 * V[j] + ((I[j] + I_A[j]) * Ra * 1e-9));

        //----- NEP ----------------------------------------------

        dPT = std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta))
              - std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta));

        dPdT = std::pow(Vg, 2) / Rn * 1e12 * dPT / (2.0 * dT); //pW/K

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
               * (std::get<0>(PowerCoolInt(DeltaT, V[j + 1] / Vg, tau, tauE))
                  - std::get<0>(PowerCoolInt(DeltaT, V[j - 1] / Vg, tau, tauE)))
               / (2 * dV); //pW/V

        G_NIS = dPdT;
        G_e = 5 * Z * Vol * std::pow(Te, 4) * 1e3; //pW/K

        G = G_e + 2 * (G_NIS - dIdT / dIdV * dPdV); //pW/K, '2' for 2 SINs

        Sv = -2 * dIdT / dIdV / G; //V/pW, for 1 bolo

        Sv /= MP;

        NEPep = 10 * EL * KBOL * Z * Vol * (std::pow(Tp, TephPOW) + std::pow(Te, TephPOW)) * 1e3
                * 1e12; //^2 pW^2/Hz

        NoiA = std::pow(vn, 2) + std::pow(in * (2 * 1e9 / dIdV + Ra) * M / MP, 2); //(V/sqrt(Hz))^2

        NEPa = NoiA / std::pow(Sv, 2); //pW^2/Hz

        //----- NEP SIN approximation ----------------------------

        dI = 2 * EL * std::abs(I[j]) / std::pow(dIdV * Sv, 2) * 1e9; //pW^2/Hz

        dPdI = 2 * 2 * EL * Pcool / (dIdV * Sv)
               * 1e9; //pW^2/Hz, second '2' is from comparison with integral

        mm = std::log(std::sqrt(2 * M_PI * KBOL * Te * Vg) / (2 * std::abs(I[j]) * Rn * 1e-9));

        dP = (0.5 + mm * mm) * (KBOL * Te) * (KBOL * Te) * std::abs(I[j]) * EL * 1e-9
             * 1e24; //pW^2/Hz

        NEPs = 2 * (dI - 2 * dPdI + dP); //pW^2/Hz, '2' for 2 SINs, all terms positive

        //fprintf(file_G, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);

        //----- NEP SIN integral ---------------------------------
        /* 
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / std::pow(dIdV * Sv, 2) * 1e9;			//pW^2/Hz
  
		dPdI = dPdI / (dIdV * Sv) * EL * std::pow(Vg, 2) / Rn * 1e12 * 1e9;	//pW^2/Hz
  
		dP *= std::pow(Vg, 3) / Rn * EL * 1e24;			//pW^2/Hz 
        
		NEPs = 2 * (dI - 2 * dPdI + dP);				//pW^2/Hz, '2' for 2 SINs, all terms positive
  
		fprintf(f9, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);
*/
        //--------------------------------------------------------

        NEPph = std::sqrt(M * MP * 2 * Pbg * 0 * 1e9 * HPLANCK * 1e-12
                          + std::pow(M * MP * Pbg * 1e-12, 2) / 1e3 / 1e9); //W/sqrt(Hz), at 0 GHz

        //NEPph = sqrt(M * MP * 2 * Pbg * 350 * 1e9 * HPLANCK * 1e-12 + pow(M * MP * Pbg * 1e-12, 2) / 1.552 / 1e9);	//W/sqrt(Hz), at 350 GHz

        NEPph *= 1e12; //pW

        NEP = std::sqrt(M * MP * (NEPep + NEPs) + NEPa + std::pow(NEPph, 2)); //all squares

        fprintf(file_Noise,
                "%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                M * (2 * V[j] + I[j] * Ra * 1e-9),
                std::sqrt(M * MP * NEPep) * 1e9 * std::abs(Sv),
                std::sqrt(M * MP * NEPs) * 1e9 * std::abs(Sv),
                std::sqrt(NoiA) * 1e9,
                NEP * 1e9 * std::abs(Sv),
                NEPph * 1e9 * std::abs(Sv),
                1e9 * std::abs(Sv) * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)));

        fprintf(file_NEP,
                "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
                M * (2 * V[j] + I[j] * Ra * 1e-9),
                MP * 1e-9 * (I[j]),
                std::sqrt(M * MP * NEPep) * 1e-12,
                std::sqrt(M * MP * NEPs) * 1e-12,
                std::sqrt(NEPa) * 1e-12,
                NEP * 1e-12,
                NEPph * 1e-12,
                std::abs(Sv) * 1e12,
                1e-12 * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)));

        fprintf(file_G, "%g\t%g\t%g\n", M * (2 * V[j] + I[j] * Ra * 1e-9), G_e, G_NIS);

        printf("Sv: %g Te: %g NEPs: %g NEPt: %g\n\n",
               std::abs(Sv) * 1e12,
               Te,
               std::sqrt(M * MP * NEPs) * 1e-12,
               NEP * 1e-12);

        //system("cls");
    }

    fclose(file_Noise);
    fclose(file_Te);
    fclose(file_NEP);
    fclose(file_G);

    double sTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    printf("\nTime spent: %f sec.", sTime);

    return NV - 1;
}
