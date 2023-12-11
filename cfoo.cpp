#include <filesystem>
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

CFoo::CFoo(size_t parnum) {
    /*
     *  Set the base parameters for fitting
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
            std::clog << parname << " = " << par[i] << ',' << SEP
                    << "to fit = " << std::boolalpha << ToFit[i]
                    << std::endl;
        }
        parfile.close();
    } catch (std::exception& e) {
        std::cerr << "Error while reading parameters: " << e.what() << std::endl;
        getchar();
        exit(0);
    }

    GetExp(fname, Iexp, Vexp, false);
    size_t countnum = CEB_2eq_parallel_lite();
    Irex.resize(countnum);
    Vrex.resize(countnum);

    Resample(Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    writeConverg(par[iParNum], ChiSq(Inum, Irex), time(nullptr));
}

double CFoo::operator()(const double dParam) {
    par[iParNum] = dParam;

    CEB_2eq_parallel_lite();
    Resample(Irex, Vrex, Iexp, Vexp, Inum, Vnum);

    return ChiSqDer(Vnum, Inum, Irex);
}

void CFoo::SeqFit(size_t iRunCount) {
    /*
     *  Fit using Golden method
     */

    Golden method(1e-3);

    std::random_device r;
    std::default_random_engine generator(r());

    for (size_t i = 0; i < iRunCount; i++) {
        std::vector<int> ParSeq(iNumParams);
        std::vector<int> Tmp(iNumParams);

        // fill `Tmp` with a sequence of ints from 0 to (iNumParams - 1)
        std::iota(Tmp.begin(), Tmp.end(), 0);

        for (size_t j = 0; j < iNumParams; ++j) {
            std::uniform_int_distribution<size_t> distribution(0, iNumParams - j);
            const size_t iRandom = distribution(generator);
            ParSeq[j] = Tmp[iRandom];
            // shift items of `Tmp` from `iRandom` to the end
            for (size_t k = iRandom + 1; k < iNumParams; ++k) {
                Tmp[k - 1] = Tmp[k];
            }
        }

        for (int j = 0; j < iNumParams; ++j)
            if (ToFit[ParSeq[j]]) {
                iParNum = ParSeq[j];
                method.ax = 0.5 * par[iParNum];
                method.bx = 1.0 * par[iParNum];
                method.cx = 2.0 * par[iParNum];
                //printf("Bracketed at: %lf %lf %lf", brent.ax, brent.bx, brent.cx);
                //scanf("%lf", &tmp);

                par[iParNum] = method.minimize(*this);
            }
        if (std::fstream params("fitparameters_new.txt", std::fstream::app); params) {
            for (const auto p: par)
                params << p << SEP;
            params << method.fmin << std::endl;
            params.close();
        }
    }
}

size_t CFoo::CEB_2eq_parallel_lite() {
    double I0; // units of current

    double Vstr; // initial voltage
    double Vfin; // final voltage

    double tau, Vg, tauE, T1, T2; // dimentionless variables

    double Pe_p, Pabs, Pleak, Pheat, Pcool, P, Ps, Pand; // for power

    double G, dPdT, dIdT, dIdV, dPdV, Sv, dPT, mm; // for noise

    double dP, dI, dPdI; // for power integrals

    double NEP, NEPs, NEPep, NEPa, NEPph; // noise equivalent power

    double DeltaT;

    double NoiA; // amplifier noise

    double G_NIS, G_e; // thermal conductivity

    time_t start = time(nullptr);

    // --------- known/guessed physical parameters

    // `bolometersInSeries` and `bolometersInParallel` may be of an integer type,
    // but as they're used in floating-point operations, let them be `double`
    auto bolometersInSeries = par[15]; // number of bolometers in series
    auto bolometersInParallel = par[16]; // number of bolometers in parallel

    double Pbg = par[0]; // incoming power for all structure [pW]
    double beta = par[1]; // returning power ratio, <1
    double TephPOW = par[2]; // exponent for Te-ph
    double Vol = par[4]; // volume of the absorber [um³]
    double Z = par[6]; // heat exchange in normal metal [nW/(K⁵×um³)]
    double Tc = par[8]; // critical temperature [K]
    double Rn = par[9] * bolometersInParallel / bolometersInSeries; // normal resistance for 1 bolometer [Ohm]
    double Rleak = par[10] * bolometersInParallel / bolometersInSeries; // leakage resistance per 1 bolometer [Ohm]
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

    double dPbg = Pbg / (bolometersInSeries * bolometersInParallel); // incoming power per 1 bolometer [pW]

    double Delta = 1.764 * Tc; // energy gap [K], Vg[eV] = Tc * 1.764 * 86.25e-6

    // if there is a file named “Te.txt”, backup its content into “Te_old.txt”
    if (std::filesystem::exists("Te.txt")) {
        std::filesystem::rename("Te.txt", "Te_old.txt");
    }

    //---------- normalized constants

    Rn = (Rn - Ra) / 2.0;
    // Rleak = Rn / gamma;							//Ohm, is not needed if Gamma<>0

    I0 = 1e9 * (Delta / Rn * K); // nA

    Vg = Delta * K;

    std::clog << "Vg = " << Vg << std::endl;

    tau = Ts / Delta;
    tauE = Te / Delta;

    //---------- calculation parameters

    Vfin = dVFinVg * Vg;
    Vstr = dVStartVg * Vg;

    auto NV = static_cast<size_t>(std::round(Vfin / dV)); // the numbers of voltage steps

    Inum.resize(NV - 1);
    Vnum.resize(NV - 1);

    std::valarray<double> V(NV + 1); // [V]
    std::iota(begin(V), end(V), 0.0);
    V = Vstr + (V * dV);

    std::ofstream file_Noise("Noise.dat");
    std::ofstream file_Te("Te.txt");
    std::ofstream file_NEP("NEP.dat");
    std::ofstream file_G("G.txt");

    file_Noise
            << "Voltage" << SEP
            << "NOISEep" << SEP
            << "NOISEs" << SEP
            << "NOISEa" << SEP
            << "NOISE" << SEP
            << "NOISEph" << SEP
            << "NOISE^2-NOISEph^2" << std::endl;
    file_Te
            << "Voltage" << SEP
            << "Current" << SEP
            << "Iqp" << SEP
            << "Iand" << SEP
            << "V/Rleak" << SEP
            << "Te" << SEP
            << "Ts" << SEP
            << "DeltaT" << SEP
            << "Pand" << SEP
            << "Pleak" << SEP
            << "Pabs" << SEP
            << "Pcool" << std::endl;
    file_NEP
            << "Voltage" << SEP
            << "Current" << SEP
            << "NEPep" << SEP
            << "NEPs" << SEP
            << "NEPa" << SEP
            << "NEP" << SEP
            << "NEPph" << SEP
            << "Sv" << SEP
            << "NEP^2-NEPph^2";
    file_G
            << "Voltage" << SEP
            << "Ge" << SEP
            << "Gnis" << std::endl;

    std::valarray<double> I(NV + 1);
    std::valarray<double> I_A(NV + 1);

    for (size_t j = 1; j < NV; ++j) // next voltage; exclude edges
    {
        constexpr double dT = 0.005;
        for (size_t n = 0; n < 5; ++n) // next interation
        {
            T1 = 0;
            T2 = 3.0 / 1.764;
            tauE = (T1 + T2) / 2;
            DeltaT = std::sqrt(1 - pow(Ts / Tc, 3.2));

            for (size_t l = 0; l < 15; ++l) // next iteration
            {
                tauE = (T1 + T2) / 2.0;

                I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9; // [nA]

                I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0; // [nA]

                Pe_p = Z * Vol
                       * (std::pow(Tp, TephPOW /*7,6,5*/)
                          - std::pow(tauE * Delta, TephPOW /*7,6,5*/))
                       * 1e3; // [pW]

                Pabs = std::pow(I[j], 2) * Ra * 1e-6; // [pW]

                Pleak = 2.0 * std::pow(V[j], 2) / Rleak * 1e12; // [pW], 2 for 2 SINs

                Pand = std::pow(I_A[j] * 1e-3 /*[uA]*/, 2) * Ra /*[Ohm]*/
                       + 2.0 * (I_A[j] * 1e3) /*[pA]*/ * V[j] /*[V]*/; // [pW], absorber + Andreev

                std::tie(Pcool, Ps) = PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE);

                Pcool *= std::pow(Vg, 2) / Rn * 1e12; // pW
                Ps *= std::pow(Vg, 2) / Rn * 1e12; // returning power from S to N

                Pheat = Pe_p + Pabs + Pand + dPbg + 2.0 * beta * Ps + Pleak;

                P = Pheat - 2.0 * Pcool;

                if (P < 0.0) {
                    T2 = tauE;
                } else {
                    T1 = tauE;
                }
            }
        }

        Te = tauE * Delta;
        Ts = tau * Delta;

        DeltaT = std::sqrt(1.0 - std::pow(Ts / Tc, 3.2));

        I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + V[j] / Rleak * 1e9;

        I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0;

        file_Te
                << (2.0 * V[j] + 1e-9 * (I[j] + I_A[j]) * Ra) * bolometersInSeries << SEP
                << 1e-9 * (I[j] + I_A[j]) * bolometersInParallel << SEP
                << 1e-9 * I[j] * bolometersInParallel << SEP
                << 1e-9 * I_A[j] * bolometersInParallel << SEP
                << 1e9 * (V[j] / Rleak) * bolometersInParallel << SEP
                << Te << SEP
                << Ts << SEP
                << DeltaT << SEP
                << Pand << SEP
                << Pleak << SEP
                << Pabs << SEP
                << Pcool << std::endl;

        std::clog << "Voltage: " << (2.0 * V[j] + 1e-9 * (I[j] + I_A[j]) * Ra) * bolometersInSeries << SEP
                << "Current: " << 1e-9 * (I[j] + I_A[j]) * bolometersInParallel << std::endl;

        //fprintf(file_NEP,"%f %f\n", M * (2 * V[j] + I[j] * Ra * 1e-9), 2 * Rsg);

        Inum[j - 1] = 1e-9 * (I[j] + I_A[j]) * bolometersInParallel;
        Vnum[j - 1] = (2.0 * V[j] + 1e-9 * (I[j] + I_A[j]) * Ra) * bolometersInSeries;

        //----- NEP ----------------------------------------------

        dPT = std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta))
              - std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta));

        dPdT = 1e12 * (std::pow(Vg, 2) / Rn) * dPT / (2.0 * dT); // [pW/K]

        dIdT = I0
               * (currentInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta)
                  - currentInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta))
               / (2.0 * dT); // [nA/K]

        dIdV = (I0
                * (currentInt(DeltaT, V[j + 1] / Vg, tau, tauE)
                   + ii * AndCurrent(DeltaT, V[j + 1] / Vg, tauE, Wt, tm)
                   - currentInt(DeltaT, V[j - 1] / Vg, tau, tauE)
                   - ii * AndCurrent(DeltaT, V[j - 1] / Vg, tauE, Wt, tm)))
               / (2.0 * dV); // [nA/V]

        dPdV = Vg * Vg / Rn * 1e12
               * (std::get<0>(PowerCoolInt(DeltaT, V[j + 1] / Vg, tau, tauE))
                  - std::get<0>(PowerCoolInt(DeltaT, V[j - 1] / Vg, tau, tauE)))
               / (2.0 * dV); // [pW/V]

        G_NIS = dPdT;
        G_e = 5.0 * Z * Vol * std::pow(Te, 4) * 1e3; // [pW/K]

        G = G_e + 2.0 * (G_NIS - dIdT / dIdV * dPdV); // [pW/K], '2' for 2 SINs

        Sv = -2.0 * dIdT / dIdV / G; // [V/pW], for 1 bolo

        Sv /= bolometersInParallel;

        NEPep = 10.0 * E * K * Z * Vol * (std::pow(Tp, TephPOW) + std::pow(Te, TephPOW)) * 1e3
                * 1e12; //^2 [pW²/Hz]

        NoiA = std::pow(vn, 2)
               + std::pow(in * (2.0 * 1e9 / dIdV + Ra) * bolometersInSeries / bolometersInParallel, 2); // [V²/Hz]

        NEPa = NoiA / std::pow(Sv, 2); // [pW²/Hz]

        //----- NEP SIN approximation ----------------------------

        dI = 1e9 * (2.0 * E * std::abs(I[j]) / std::pow(dIdV * Sv, 2)); // [pW²/Hz]

        dPdI = 1e9 * (2.0 * 2.0 * E * Pcool / (dIdV * Sv)); // [pW^2/Hz], second '2' is from comparison with integral

        mm = std::log(std::sqrt(2.0 * M_PI * K * Te * Vg) / (2.0 * std::abs(I[j]) * Rn * 1e-9));

        dP = (0.5 + mm * mm) * (K * Te) * (K * Te) * std::abs(I[j]) * E * 1e-9
             * 1e24; // [pW²/Hz]

        NEPs = 2.0 * (dI - 2.0 * dPdI + dP); // [pW²/Hz], '2' for 2 SINs, all terms positive

        //fprintf(file_G, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);

        //----- NEP SIN integral ---------------------------------
        /*
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / std::pow(dIdV * Sv, 2) * 1e9;			// [pW²/Hz]

		dPdI = dPdI / (dIdV * Sv) * EL * std::pow(Vg, 2) / Rn * 1e12 * 1e9;	// [pW²/Hz]

		dP *= std::pow(Vg, 3) / Rn * EL * 1e24;			// [pW²/Hz]

		NEPs = 2 * (dI - 2 * dPdI + dP);				// [pW²/Hz], '2' for 2 SINs, all terms positive

		fprintf(f9, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);
        */
        //--------------------------------------------------------

        NEPph = std::sqrt(std::pow(1e-12 * Pbg * bolometersInSeries * bolometersInParallel, 2) / 1e3 / 1e9);
        // [W/sqrt(Hz)], at 0 GHz

        //NEPph = sqrt(M * MP * 2 * Pbg * 350 * 1e9 * HPLANCK * 1e-12 + pow(M * MP * Pbg * 1e-12, 2) / 1.552 / 1e9);	//W/sqrt(Hz), at 350 GHz

        NEPph *= 1e12; // [pW]

        NEP = std::sqrt((NEPep + NEPs) * bolometersInSeries * bolometersInParallel + NEPa + std::pow(NEPph, 2));
        //all squares

        file_Noise
                << (2.0 * V[j] + 1e-9 * I[j] * Ra) * bolometersInSeries << SEP
                << 1e9 * std::sqrt(NEPep * (bolometersInSeries * bolometersInParallel)) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NEPs * (bolometersInSeries * bolometersInParallel)) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NoiA) << SEP
                << 1e9 * NEP * std::abs(Sv) << SEP
                << 1e9 * NEPph * std::abs(Sv) << SEP
                << 1e9 * std::abs(Sv) * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)) << std::endl;

        file_NEP
                << (2.0 * V[j] + 1e-9 * I[j] * Ra) * bolometersInSeries << SEP
                << 1e-9 * I[j] * bolometersInParallel << SEP
                << 1e-12 * std::sqrt(NEPep * (bolometersInSeries * bolometersInParallel)) << SEP
                << 1e-12 * std::sqrt(NEPs * (bolometersInSeries * bolometersInParallel)) << SEP
                << 1e-12 * std::sqrt(NEPa) << SEP
                << 1e-12 * NEP << SEP
                << 1e-12 * NEPph << SEP
                << 1e12 * std::abs(Sv) << SEP
                << 1e-12 * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)) << std::endl;

        file_G
                << (2.0 * V[j] + 1e-9 * (I[j] * Ra)) * bolometersInSeries << SEP
                << G_e << SEP
                << G_NIS << std::endl;

        std::clog
                << "Sv: " << 1e12 * std::abs(Sv) << SEP
                << "Te: " << Te << SEP
                << "NEPs: " << 1e-12 * std::sqrt(NEPs * (bolometersInSeries * bolometersInParallel)) << SEP
                << "NEPt: " << 1e-12 * NEP << SEP
                << '\n' << std::endl;
    }

    file_Noise.close();
    file_Te.close();
    file_NEP.close();
    file_G.close();

    std::clog << '\n' << "Time spent: " << difftime(time(nullptr), start) << " sec." << std::endl;

    return NV - 1;
}
