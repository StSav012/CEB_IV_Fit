// Defines the entry point for the console application.
//Pleak is included in Pcool

//for bolometers (absorber + 2SINs) w/ Andreev Current
#ifndef CEB_2EQ_PARALLEL_LITE_H
#define CEB_2EQ_PARALLEL_LITE_H

#include <cmath>
#include <ctime>
#include <tuple>

#include "constants.h"

/*
 * Never used
 */
inline double current(double v, double tau) {
    /*  Compute the current using approximation [delta(0) / eRn]
     *
     *  `v` is voltage,
     *  `tau` is for temperature
     */

    double s;

    const double a0 = 1 + 0.375 * tau - 0.1171875 * tau * tau;
    const double a1 = tau * a0 * a0;
    double a2 = 1 + std::exp((std::abs(v) - 1) / tau - (1.15 + tau));
    double a3 = (v * v - 1) / (1 + std::exp(-(v * v - 1) / tau));

    if (v > 0)
        s = 1;
    else
        s = -1;

    // [delta(0) / eRn]
    return s * std::sqrt(2.0 * M_PI * a1 / a2 + a3)
           * (1.0 / (2.0 * std::exp((1.0 - v) / tau) + 1.0)
              + 1.0 / (2.0 * std::exp((1.0 + v) / tau) + 1.0));
}

inline double currentInt(const double DT, const double v, const double tau, double const tauE) {
    /*
     * Compute the current using exact integral, rectangle method
     *
     * `DT` is delta (energy gap),
     * `v` is voltage,
     * `tau` is for temperature,
     * `tauE` is for electron temperature
    */

    const double dx = INTEGRATION_SCALE * DT;
    const double inf = ESSENTIALLY_INFINITY_SCALE * DT;
    const double DT2 = std::pow(DT, 2);
    const double v1 = std::abs(v);

    double pos_i = 0, neg_i = 0;

    double x = DT; // energy
    while (x <= inf) {
        constexpr double accuracy = 1e-6;
        const double a2 = std::sqrt(std::pow(x, 2) - DT2);

        double a0 = 1.0 / (std::exp((x - v1) / tauE) + 1);
        double a1 = 1.0 / (std::exp(x / tau) + 1);
        double i_step = x / a2 * (a0 - a1);
        // accuracy check at every iteration
        if (std::abs(i_step) < accuracy)
            break;
        pos_i += i_step;

        a0 = 1.0 / (std::exp((-x - v1) / tauE) + 1);
        a1 = 1.0 / (std::exp(-x / tau) + 1);
        i_step = x / a2 * (a0 - a1);
        // accuracy check at every iteration
        if (std::abs(i_step) < accuracy)
            break;
        neg_i += i_step;

        x += dx;
    }

    double s;
    if (v > 0)
        s = 1;
    else
        s = -1;

    return (pos_i + neg_i) * dx * s;
}

inline double PowerCoolPiece(const double v, const double tau) {
    /*
     *  A part of `PowerCool`
     */
    const double a0 = (1 - v) / tau;
    const double tmp1 = 2.0 * std::exp(a0);
    const double tmp2 = std::exp(2.5 * (a0 + 2));
    const double PItau = M_PI * tau;
    const double a1 = std::sqrt(2.0 * PItau) * ((1 - v) / (tmp1 + 1.28) + 0.5 * tau / (tmp1 + 0.64))
                      / (1.0 / tmp2 + 1.0);
    double a3;
    if (v - 1 - tau > 0) {
        const double a2 = std::sqrt(std::pow(v, 2) - 1.0);
        a3 = 0.5 * (-v * a2 + std::log(v + a2) + std::pow(PItau, 2) / 3.0 * v / a2) / (tmp2 + 1.0);
    } else {
        a3 = 0.0;
    }

    return a1 + a3;
}

inline double PowerCool(const double v, const double tau, const double tauE) {
    /*
     *  Compute the cooling power of a SIN junction
     *
     *  `v` is voltage,
     *  `tau` is for temperature,
     *  `tauE` is for electron temperature
     */

    const double a1 = PowerCoolPiece(v, tauE);
    const double b1 = PowerCoolPiece(-v, tauE);
    const double c1 = PowerCoolPiece(0.0, tau);

    return a1 + b1 - c1;
}

inline double AndCurrent(
    const double DT, const double v, const double tauE, const double Wt, const double tm) {
    /*
     *  Compute the Andreev current using exact integral, rectangle method
     *
     *  `DT` is delta (energy gap),
     *  `v` is voltage,
     *  `tauE` is for electron temperature,
     *  `Wt` is omega with ~ from Vasenko 2010, for the transparency of the barrier,
     *  `tm` is for the depairing energy
     */

    const double dx = INTEGRATION_SCALE * DT;
    const double DT2 = std::pow(DT, 2);
    const double twoWt = 2.0 * Wt;
    const double twoTauE = 2.0 * tauE;

    double i = 0.0;

    //FILE *f9=fopen("SINi.dat","w");

    double x = 0.0; //energy
    while (x <= DT) {
        const double DT2x2 = DT2 - std::pow(x, 2);
        const double sqrtDT2x2 = std::sqrt(DT2x2);
        const double da0a1 = std::tanh((x + v) / twoTauE) - std::tanh((x - v) / twoTauE);
        const double a2 = twoWt * sqrtDT2x2 / tm
                          / (std::pow(x * (twoWt - sqrtDT2x2 / DT), 2) + DT2x2 / std::pow(DT * tm, 2));
        i += DT / sqrtDT2x2 * da0a1 * a2;
        //fprintf(f9,"%g %g\n", x, i);
        x += dx;
    }
    i *= dx;
    return i;
    //fclose(f9);
}

inline auto PowerCoolInt(const double DT, const double v, const double tau, const double tauE) {
    /*
     *  Compute the integral of the cooling power of SIN junction, rectangle method
     *
     *  `DT` is delta (energy gap),
     *  `v` is voltage,
     *  `tau` is for temperature,
     *  `tauE` is for electron temperature
     */

    //E = x, V = v, T = tau, [x] = [v] = [tau]

    const double dx = INTEGRATION_SCALE * DT;
    const double inf = ESSENTIALLY_INFINITY_SCALE * DT;
    const double v1 = std::abs(v);
    const double DT2 = std::pow(DT, 2);

    double po = 0.0; // power
    double ps = 0.0; // SIN cooling power

    //FILE *f8 = fopen("DOS.txt","w");  //density of states

    double x = DT; // energy
    while (x <= inf) {
        double x2 = std::pow(x, 2);
        double a = std::sqrt(x2 - DT2);

        double da2a1 = (1.0 / (std::exp((x - v1) / tauE) + 1.0) - 1.0 / (std::exp(x / tau) + 1.0)) / a;
        po += x * (x - v1) * da2a1;
        ps += x2 * da2a1;

        da2a1 = (1.0 / (std::exp((-x - v1) / tauE) + 1.0) - 1.0 / (std::exp(-x / tau) + 1.0)) / a;
        po += x * (-x - v1) * da2a1;
        ps += -x2 * da2a1;

        x += dx;
    }
    po *= dx;
    ps *= dx;
    return std::make_tuple(po, ps);
}
#endif // CEB_2EQ_PARALLEL_LITE_H
