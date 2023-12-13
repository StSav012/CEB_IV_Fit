#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "constants.h"
#include "utils.h"
#include "CEBNumericModel.h"
#include "MinimizationAlgorithms.h"

#include "IVParamFitter.h"

IVParamFitter::IVParamFitter() {
    /*
     *  Set the base parameters for fitting
     */
    if (std::ifstream parfile("startparams.txt"); parfile) {
        std::string parname;
        double parvalue;
        bool parfit;
        while (!parfile.eof()) {
            parfile >> std::skipws >> parname >> parvalue >> parfit;

            par[parname] = parvalue;
            ToFit[parname] = parfit;

            std::clog << std::left << std::setw(18) << std::format("{} = {},", parname, parvalue)
                    << std::internal << "to fit = " << std::boolalpha << parfit
                    << std::endl;
        }
        parfile.close();
    } else {
        throw std::runtime_error("Can't read \"startparams.txt\"");
    }
}

double IVParamFitter::operator()(const double dParam) {
    par[parameterName] = dParam;

    CEB_2eq_parallel_lite();
    auto [Irex, Vrex] = resample();

    return ChiSqDer(Vnum, Inum, Irex);
}

void IVParamFitter::SeqFit(const size_t runCount, const std::valarray<double>& Irex) {
    /*
     *  Fit using Golden method
     */

    std::random_device r;
    std::default_random_engine generator(r());

    writeConverg(par[parameterName], ChiSq(Inum, Irex), std::chrono::steady_clock::now());

    for (size_t run = 0; run < runCount; ++run) {
        std::clog << "SeqFit run " << run << std::endl;

        std::vector<std::string> ParSeq;
        for (const auto& [key, value]: ToFit) {
            if (value) {
                ParSeq.push_back(key);
            }
        }
        std::ranges::shuffle(ParSeq, generator);

        double fmin = NAN;
        for (auto& parname: ParSeq) {
            parameterName = parname;
            std::tie(par[parameterName], fmin) = GoldenMinimize(
                *this,
                0.5 * par[parameterName], 2.0 * par[parameterName],
                1.0 * par[parameterName],
                1e-3
            );
        }

        // store the parameters after minimization and the result
        bool appendNewLine = std::filesystem::exists("fitparameters_new.txt")
                             && std::filesystem::file_size("fitparameters_new.txt");
        if (std::fstream params("fitparameters_new.txt", std::ios::app); params) {
            if (appendNewLine) {
                params << std::endl;
            }
            params << std::format("time = {}", std::chrono::system_clock::now()) << std::endl;
            for (const auto& [parname, parvalue]: par)
                params << std::format("{} = {} ({})", parname, parvalue,
                                      ToFit[parname] ? std::string("fit") : std::string("skip")) << std::endl;
            params << std::format("fmin = {}", fmin) << std::endl;
            params.close();
        } else {
            throw std::runtime_error("Unable to append to \"fitparameters_new.txt\"");
        }
    }
}

size_t IVParamFitter::loadExperimentData(const std::string& filename, const bool removeOffset) {
    std::tie(Iexp, Vexp) = getExperimentalData(filename, removeOffset);
    if (Iexp.size() != Vexp.size()) {
        throw std::length_error("Experimental I and V must be of the same size");
    }
    return Iexp.size();
}

std::tuple<std::valarray<double>, std::valarray<double>> IVParamFitter::resample() const {
    return Resample(Iexp, Vexp, Inum, Vnum);
}

size_t IVParamFitter::CEB_2eq_parallel_lite() {
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

    // --------- known/guessed physical parameters

    // `bolometersInSeries` and `bolometersInParallel` may be of an integer type,
    // but as they're used in floating-point operations, let them be `double`
    const auto bolometersInSeries = par["M"]; // number of bolometers in series
    const auto bolometersInParallel = par["MP"]; // number of bolometers in parallel
    const auto totalBolometersNumber = bolometersInSeries * bolometersInParallel;

    // incoming power for all structure [pW]
    const double Pbg = par["Pbg"];
    // returning power ratio, <1
    const double beta = par["beta"];
    // exponent for Te-ph, 7, 6, or 5
    const double TephPOW = par["TephPOW"];
    // volume of the absorber [um³]
    const double Vol = par["Vol"];
    // heat exchange in normal metal [nW/(K⁵×um³)]
    const double Z = par["Z"];
    // critical temperature [K]
    const double Tc = par["Tc"];
    // normal resistance for 1 bolometer [Ohm]
    double Rn = par["Rn"] * bolometersInParallel / bolometersInSeries;
    // leakage resistance per 1 bolometer [Ohm]
    const double Rleak = par["Rleak"] * bolometersInParallel / bolometersInSeries;
    // transparency of the barrier
    const double Wt = par["Wt"];
    // depairing energy
    const double tm = par["tm"];
    // coefficient for Andreev current
    const double ii = par["ii"];
    // normal resistance of 1 absorber [Ohm]
    const double Ra = par["Ra"];
    // phonon temperature [K]
    const double Tp = par["Tp"];
    // voltage range end [V]
    const double dVFinVg = par["dVFinVg"];
    // voltage range start [V]
    const double dVStartVg = par["dVStartVg"];
    // voltage step [V]
    const double dV = par["dV"];

    // electron temperature to be found [K]
    double Te = Tp;
    // electron temperature in superconductor [K]
    const double Ts = Tp;

    const double DeltaT = std::sqrt(1.0 - std::pow(Ts / Tc, 3.2));
    // incoming power per 1 bolometer [pW]
    const double dPbg = Pbg / totalBolometersNumber;
    // energy gap [K], Vg[eV] = Tc * 1.764 * 86.25e-6
    const double Delta = 1.764 * Tc;

    // if there is a file named “Te.txt”, backup its content into “Te_old.txt”
    if (std::filesystem::exists("Te.txt")) {
        std::filesystem::rename("Te.txt", "Te_old.txt");
    }

    //---------- normalized constants

    Rn = (Rn - Ra) / 2.0;

    const double I0 = 1e9 * (Delta / Rn * K); // [nA], units of current

    const double Vg = Delta * K; // dimentionless
    std::clog << "Vg = " << Vg << std::endl;

    const double tau = Ts / Delta; // dimentionless
    double tauE = Te / Delta; // dimentionless

    //---------- calculation parameters

    // initial voltage
    const double Vfin = dVFinVg * Vg;
    // final voltage
    const double Vstr = dVStartVg * Vg;

    const auto NV = static_cast<size_t>(std::round((Vfin - Vstr) / dV)); // the number of voltage steps
    if (!NV) {
        throw std::length_error("No voltage steps to do");
    }

    Inum.resize(NV - 1);
    Vnum.resize(NV - 1);

    std::valarray<double> V(NV + 1); // [V]
    std::iota(begin(V), end(V), 0.0);
    V = Vstr + (V * dV);

    std::ofstream file_Noise("Noise.txt");
    if (!file_Noise) {
        throw std::runtime_error("Unable to write \"file_Noise.txt\"");
    }
    std::ofstream file_Te("Te.txt");
    if (!file_Te) {
        throw std::runtime_error("Unable to write \"file_Te.txt\"");
    }
    std::ofstream file_NEP("NEP.txt");
    if (!file_NEP) {
        throw std::runtime_error("Unable to write \"NEP.txt\"");
    }
    std::ofstream file_G("G.txt");
    if (!file_G) {
        throw std::runtime_error("Unable to write \"G.txt\"");
    }

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
            << "NEP^2-NEPph^2" << std::endl;
    file_G
            << "Voltage" << SEP
            << "Ge" << SEP
            << "Gnis" << std::endl;

    std::valarray<double> I(NV + 1);
    std::valarray<double> I_A(NV + 1);

    for (size_t j = 1; j < NV; ++j) // next voltage; exclude edges
    {
        constexpr double dT = 0.005; // temperature step for derivative calculations

        double Pabs, Pleak, Pcool, Ps, Pand; // for power

        for (size_t n = 0; n < 5; ++n) // next interation
        {
            double T1 = 0.0; // dimentionless
            double T2 = 3.0 / 1.764; // dimentionless

            for (size_t l = 0; l < 15; ++l) // next iteration
            {
                tauE = (T1 + T2) / 2.0;

                I[j] = currentInt(DeltaT, V[j] / Vg, tau, tauE) * I0 + 1e9 * (V[j] / Rleak); // [nA]

                I_A[j] = ii * AndCurrent(DeltaT, V[j] / Vg, tauE, Wt, tm) * I0; // [nA]

                const double Pe_p = Z * Vol
                                    * (std::pow(Tp, TephPOW)
                                       - std::pow(tauE * Delta, TephPOW))
                                    * 1e3; // [pW]

                Pabs = std::pow(I[j], 2) * Ra * 1e-6; // [pW]

                Pleak = numberOfSINs * std::pow(V[j], 2) / Rleak * 1e12; // [pW]

                Pand = std::pow(I_A[j] * 1e-3 /*[uA]*/, 2) * Ra /*[Ohm]*/
                       + 2.0 * (I_A[j] * 1e3) /*[pA]*/ * V[j] /*[V]*/; // [pW], absorber + Andreev

                std::tie(Pcool, Ps) = PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE);

                Pcool *= std::pow(Vg, 2) / Rn * 1e12; // [pW]
                Ps *= std::pow(Vg, 2) / Rn * 1e12; // returning power from S to N

                if (const double Pheat = Pe_p + Pabs + Pand + dPbg + 2.0 * beta * Ps + Pleak; Pheat < 2.0 * Pcool) {
                    T2 = tauE;
                } else {
                    T1 = tauE;
                }
            }
        }

        Te = tauE * Delta;

        Inum[j - 1] = 1e-9 * (I[j] + I_A[j]) * bolometersInParallel;
        Vnum[j - 1] = (2.0 * V[j] + 1e-9 * (I[j] + I_A[j]) * Ra) * bolometersInSeries;

        file_Te
                << Vnum[j - 1] << SEP
                << Inum[j - 1] << SEP
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

        //----- NEP ----------------------------------------------

        const double dPT = std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta))
                           - std::get<0>(PowerCoolInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta));

        const double dPdT = 1e12 * (std::pow(Vg, 2) / Rn) * dPT / (2.0 * dT); // [pW/K]

        const double dIdT = I0
                            * (currentInt(DeltaT, V[j] / Vg, tau, tauE + dT / Delta)
                               - currentInt(DeltaT, V[j] / Vg, tau, tauE - dT / Delta))
                            / (2.0 * dT); // [nA/K]

        const double dIdV = I0
                            * (currentInt(DeltaT, V[j + 1] / Vg, tau, tauE)
                               + ii * AndCurrent(DeltaT, V[j + 1] / Vg, tauE, Wt, tm)
                               - currentInt(DeltaT, V[j - 1] / Vg, tau, tauE)
                               - ii * AndCurrent(DeltaT, V[j - 1] / Vg, tauE, Wt, tm))
                            / (2.0 * dV); // [nA/V]

        const double dPdV = std::pow(Vg, 2) / Rn * 1e12
                            * (std::get<0>(PowerCoolInt(DeltaT, V[j + 1] / Vg, tau, tauE))
                               - std::get<0>(PowerCoolInt(DeltaT, V[j - 1] / Vg, tau, tauE)))
                            / (2.0 * dV); // [pW/V]

        // thermal conductivity
        const double G_NIS = dPdT;
        const double G_e = 5.0 * Z * Vol * std::pow(Te, 4) * 1e3; // [pW/K]

        const double G = G_e + numberOfSINs * (G_NIS - dIdT / dIdV * dPdV); // [pW/K]

        const double Sv = -2.0 * dIdT / dIdV / G / bolometersInParallel; // [V/pW], for 1 bolo

        const double NEPep = 10.0 * E * K * Z * Vol * (std::pow(Tp, TephPOW) + std::pow(Te, TephPOW)) * 1e3
                             * 1e12; //^2 [pW²/Hz]

        // amplifier noise
        const double NoiA = std::pow(VOLTAGE_NOISE_2_AMPS, 2)
                            + std::pow(
                                CURRENT_NOISE_2_AMPS * (2.0 * 1e9 / dIdV + Ra) * bolometersInSeries /
                                bolometersInParallel,
                                2); // [V²/Hz]

        const double NEPa = NoiA / std::pow(Sv, 2); // [pW²/Hz]

        //----- NEP SIN approximation ----------------------------

        const double dI = 1e9 * (2.0 * E * std::abs(I[j]) / std::pow(dIdV * Sv, 2)); // [pW²/Hz]

        const double dPdI = 1e9 * (2.0 * 2.0 * E * Pcool / (dIdV * Sv));
        // [pW^2/Hz], second '2' is from comparison with integral

        const double mm = std::log(std::sqrt(2.0 * M_PI * K * Te * Vg) / (2.0 * std::abs(I[j]) * Rn * 1e-9));

        const double dP = (0.5 + std::pow(mm, 2)) * std::pow(K * Te, 2) * std::abs(I[j]) * E * 1e-9
                          * 1e24; // [pW²/Hz]

        const double NEPs = numberOfSINs * (dI - 2.0 * dPdI + dP); // [pW²/Hz], all terms positive

        //fprintf(file_G, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);

        //----- NEP SIN integral ---------------------------------
        /*
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / std::pow(dIdV * Sv, 2) * 1e9;			// [pW²/Hz]

		dPdI = dPdI / (dIdV * Sv) * EL * std::pow(Vg, 2) / Rn * 1e12 * 1e9;	// [pW²/Hz]

		dP *= std::pow(Vg, 3) / Rn * EL * 1e24;			// [pW²/Hz]

		NEPs = numberOfSINs * (dI - 2 * dPdI + dP);				// [pW²/Hz], all terms positive

		fprintf(f9, "%g %g %g %g\n",  M * (2 * V[j] + I[j] * Ra * 1e-9), dI, 2 * dPdI, dP);
        */
        //--------------------------------------------------------

        const double NEPph = 1e12 * std::sqrt(std::pow(1e-12 * Pbg * totalBolometersNumber, 2) / 1e3 / 1e9);
        // [pW/sqrt(Hz)], at 0 GHz

        // const double NEPph = 1e12 * std::sqrt(totalBolometersNumber * 2.0 * 1e9 * Pbg * 350.0 * 1e-12 * H + std::pow(1e-12 * Pbg * totalBolometersNumber, 2) / 1.552 / 1e9);	// [pW/sqrt(Hz)], at 350 GHz

        const double NEP = std::sqrt((NEPep + NEPs) * totalBolometersNumber + NEPa + std::pow(NEPph, 2));
        //all squares

        file_Noise
                << (2.0 * V[j] + 1e-9 * I[j] * Ra) * bolometersInSeries << SEP
                << 1e9 * std::sqrt(NEPep * totalBolometersNumber) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NEPs * totalBolometersNumber) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NoiA) << SEP
                << 1e9 * NEP * std::abs(Sv) << SEP
                << 1e9 * NEPph * std::abs(Sv) << SEP
                << 1e9 * std::abs(Sv) * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)) << std::endl;

        file_NEP
                << (2.0 * V[j] + 1e-9 * I[j] * Ra) * bolometersInSeries << SEP
                << 1e-9 * I[j] * bolometersInParallel << SEP
                << 1e-12 * std::sqrt(NEPep * totalBolometersNumber) << SEP
                << 1e-12 * std::sqrt(NEPs * totalBolometersNumber) << SEP
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
                << std::setw(static_cast<int>(std::ceil(std::log10(NV)))) << j << '/' << NV - 1 << ':' << SEP
                << "Voltage: " << std::setw(12) << Vnum[j - 1] << SEP
                << "Current: " << std::setw(12) << Inum[j - 1] << SEP
                << "Sv: " << std::setw(12) << 1e12 * std::abs(Sv) << SEP
                << "Te: " << std::setw(12) << Te << SEP
                << "NEPs: " << std::setw(12) << 1e-12 * std::sqrt(NEPs * totalBolometersNumber) << SEP
                << "NEPt: " << std::setw(12) << 1e-12 * NEP << std::endl;
    }

    file_Noise.close();
    file_Te.close();
    file_NEP.close();
    file_G.close();

    std::clog
            << "Time spent: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            << std::endl;

    return NV - 1;
}
