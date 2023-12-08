// Defines the entry point for the console application.
//Pleak is included in Pcool

//for bolometers (absorber + 2SINs) w/ Andreev Current

#include <cmath>
#include <cstdio>
#include <ctime>

#include "constants.h"

double current(float v, float tau)
{
    /*
     * Compute the current using approximation [delta(0) / eRn]
     * 
     * `v` is voltage,
     * `tau` is for temperature
    */

    double i; // [delta(0) / eRn]
    double a0, a1, a2, a3;
    float s;

    a0 = 1 + 0.375 * tau - 0.1171875 * tau * tau;
    a1 = tau * a0 * a0;
    a2 = 1 + std::exp((std::abs(v) - 1) / tau - (1.15 + tau));
    a3 = (v * v - 1) / (1 + std::exp(-(v * v - 1) / tau));

    if (v > 0)
        s = 1;
    else
        s = -1;

    i = s * std::sqrt(2 * PI * a1 / a2 + a3)
        * (1 / (2 * std::exp((1 - v) / tau) + 1) + 1 / (2 * std::exp((1 + v) / tau) + 1));
    return i;
}

double currentInt(float DT, float v, float tau, float tauE)
{
    /* 
     * Compute the current using exact integral, rectangle method
     * 
     * `DT` is delta (energy gap),
     * `v` is voltage,
     * `tau` is for temperature,
     * `tauE` is for electron temperature
    */

    double i, v1;
    double a0, a1, a2;
    double x; //energy
    float dx;
    int l, inf;
    float s;
    dx = 0.2e-3;
    inf = (int) (20 / dx);
    i = 0;
    x = DT + dx;
    v1 = std::abs(v);

    double prev_i = 0;      //previous 'i' value
    double accuracy = 1e-6; //calculation accuracy

    //#pragma omp parallel for default(shared) private(x, a0, a1, a2) reduction(+:i) num_threads(THREADS)
    for (l = 0; l <= inf; l++) {
        a0 = 1 / (std::exp((x - v1) / tauE) + 1);
        a1 = 1 / (std::exp(x / tau) + 1);
        a2 = std::sqrt(x * x - DT * DT);
        i += std::abs(x) / a2 * (a0 - a1);
        x += dx;

        //accuracy check on every iteration
        if (std::abs(i - prev_i) < accuracy)
            break;
        prev_i = i;
    }
    x = -DT - dx;

    //#pragma omp parallel for default(shared) private(x, a0, a1, a2) reduction(+:i) num_threads(THREADS)
    for (l = 0; l <= inf; l++) {
        a0 = 1 / (std::exp((x - v1) / tauE) + 1);
        a1 = 1 / (std::exp(x / tau) + 1);
        a2 = std::sqrt(x * x - DT * DT);
        i += std::abs(x) / a2 * (a0 - a1);
        x -= dx;

        //accuracy check on every iteration
        if (std::abs(i - prev_i) < accuracy)
            break;
        prev_i = i;
    }
    if (v > 0)
        s = 1;
    else
        s = -1;
    i *= dx * s;
    return i;
}

double PowerCool(float v, float tau, float tauE)
{
    /*
     * Compute the cooling power of a SIN junction
     *
     * `v` is voltage, 
     * `tau` is for temperature,
     * `tauE` is for electron temperature
    */

    double p;
    double a0, a1, a2, a3;
    double b0, b1, b2, b3;
    double c0, c1;

    a3 = 0;
    b3 = 0;
    v = std::abs(v);

    a0 = (1 - v) / tauE;
    a1 = std::sqrt(2 * PI * tauE)
         * ((1 - v) / (2 * std::exp(a0) + 1.28) + 0.5 * tauE / (2 * std::exp(a0) + 0.64))
         / (std::exp(-2.5 * (a0 + 2)) + 1);
    if ((v - 1 - tauE) > 0) {
        a2 = std::sqrt(v * v - 1);
        a3 = 0.5 * (-v * a2 + std::log(v + a2) + PI * PI * tauE * tauE / 3 * v / a2)
             / (std::exp(2.5 * (a0 + 2)) + 1);
    }

    b0 = (1 + v) / tauE;
    b1 = std::sqrt(2 * PI * tauE)
         * ((1 + v) / (2 * std::exp(b0) + 1.28) + 0.5 * tauE / (2 * std::exp(b0) + 0.64))
         / (std::exp(-2.5 * (b0 + 2)) + 1);
    if ((-v - 1 - tauE) > 0) {
        b2 = sqrt(v * v - 1);
        b3 = 0.5 * (v * b2 + std::log(-v + b2) - PI * PI * tauE * tauE / 3 * v / b2)
             / (std::exp(2.5 * (b0 + 2)) + 1);
    }

    c0 = 1 / tau;
    c1 = 2 * std::sqrt(2 * PI * tau)
         * (1 / (2 * std::exp(c0) + 1.28) + 0.5 * tau / (2 * std::exp(c0) + 0.64))
         / (std::exp(-2.5 * (c0 + 2)) + 1);

    p = (a1 + a3 + b1 + b3 - c1);
    return p;
}

double AndCurrent(float DT, float v, float tauE, float Wt, float tm)
{
    /*
     * Compute the Andreev current using exact integral, rectangle method
     * 
     * `DT` is delta (energy gap),
     * `v` is voltage,
     * `tauE` is for electron temperature,
     * `Wt` is omega with ~ from Vasenko 2010, for the transparency of the barrier,
     * `tm` is for the depairing energy
    */

    double i, v1;
    double a0, a1, a2;
    double x; //energy
    float dx;
    int l, inf;
    //FILE *f9;

    dx = 0.2e-3;
    inf = (int) (DT / dx) - 1;

    i = 0;
    x = dx;

    //f9=fopen("SINi.dat","w");

    //#pragma omp parallel for default(shared) private(a0, a1, a2) reduction(+:i) num_threads(THREADS)
    for (l = 0; l <= inf - 1; l++) {
        a0 = std::tanh((x + v) / (2 * tauE));
        a1 = std::tanh((x - v) / (2 * tauE));
        a2 = 2 * Wt * std::sqrt(DT * DT - x * x) / tm
             / (std::pow(2 * Wt * x - x * std::sqrt(DT * DT - x * x) / DT, 2)
                + (DT * DT - x * x) / std::pow(DT * tm, 2));
        i += DT / std::sqrt(DT * DT - x * x) * (a0 - a1) * a2;
        //fprintf(f9,"%g %g\n", x, i);
        x += dx;
    }
    i *= dx;
    return i;
    //fclose(f9);
}

double PowerCoolInt(float DT, float v, float tau, float tauE, double *Ps)
{
    /*
     * Compute the integral of the cooling power of SIN junction, rectangle method
     * 
     * `DT` is delta (energy gap),
     * `v` is voltage,
     * `tau` is for temperature,
     * `tauE` is for electron temperature,
     * `Ps` is for SIN cooling power
    */

    //E = x, V = v, T = tau, [x] = [v] = [tau]

    double p, q;
    double a, a1, a2;
    double b, b1;
    double x, dx, de; //energy
    int l, inf;
    double po; //power
    FILE *f8;  //density of states

    po = 0;
    *Ps = 0;
    v = std::abs(v);

    dx = 2e-4;
    inf = (int) (20 / dx);

    //f8 = fopen("DOS.txt","w");

    de = DT;
    x = de + dx;

    //#pragma omp parallel for default(shared) private(x, a, a1, a2) reduction(+:po) num_threads(THREADS)
    for (l = 0; l <= inf; l++) {
        a = std::sqrt(x * x - DT * DT);

        a2 = 1 / (std::exp((x - v) / tauE) + 1);
        a1 = 1 / (std::exp(x / tau) + 1);

        po += std::abs(x) * (x - v) * (a2 - a1) / a;
        *Ps += std::abs(x) * (x) * (a2 - a1) / a;
        x += dx;
    }
    x = -de - dx;

    //#pragma omp parallel for default(shared) private(x, a, a1, a2) reduction(+:po) num_threads(THREADS)
    for (l = 0; l <= inf; l++) {
        a = std::sqrt(x * x - DT * DT);

        a2 = 1 / (std::exp((x - v) / tauE) + 1);
        a1 = 1 / (std::exp(x / tau) + 1);

        po += std::abs(x) * (x - v) * (a2 - a1) / a;
        *Ps += std::abs(x) * x * (a2 - a1) / a;
        x -= dx;
    }
    po *= dx;
    *Ps *= dx;
    return po;
}
