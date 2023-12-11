#include <cmath>
#include <stdexcept>

#include "utils.h"

#include "coff.h"

COff::COff(const std::valarray<double>& Iexp, const std::valarray<double>& Vexp) {
    if (Iexp.size() != Vexp.size()) {
        throw std::length_error("I and V must be of the same size");
    }

    Iref = Iexp;
    Vref = Vexp;
}

double COff::operator()(const double dOffset) const {
    /*
     *  Make an object of the class act like a callable
     *  
     *  TODO: describe what the call is doing.
     */
    size_t lLowI = 0;
    double tmp = std::abs(Iref[0]);
    const std::valarray<double> Vofx = Vref - dOffset;
    const std::valarray<double> Iofx = Iref;
    const size_t count = Iref.size();
    for (size_t i = 0; i < count; i++) {
        if (std::abs(Iref[i]) < tmp) {
            lLowI = i;
            tmp = std::abs(Iref[i]);
        }
    }
    if (Iofx[lLowI] < 0)
        ++lLowI;
    const size_t lLengthPos = count - lLowI;
    const size_t lLengthNeg = count - lLengthPos;
    std::valarray<double> Vlow, Ilow, Vhigh, Ihigh;
    size_t lLengthLow, lLengthHigh;
    if (lLengthPos > lLengthNeg) {
        lLengthHigh = lLengthPos;
        lLengthLow = lLengthNeg;
        Vhigh.resize(lLengthHigh);
        Ihigh.resize(lLengthHigh);
        Vlow.resize(lLengthLow);
        Ilow.resize(lLengthLow);
        for (size_t i = 0; i < lLengthHigh; i++) {
            Vhigh[i] = std::abs(Vofx[lLowI + i]);
            Ihigh[i] = std::abs(Iofx[lLowI + i]);
        }
        for (size_t i = 0; i < lLengthLow; i++) {
            Vlow[i] = std::abs(Vofx[lLowI - 1 - i]);
            Ilow[i] = std::abs(Iofx[lLowI - 1 - i]);
        }
    } else {
        lLengthHigh = lLengthNeg;
        lLengthLow = lLengthPos;
        Vhigh.resize(lLengthHigh);
        Ihigh.resize(lLengthHigh);
        Vlow.resize(lLengthLow);
        Ilow.resize(lLengthLow);
        for (size_t i = 0; i < lLengthHigh; i++) {
            Vhigh[i] = std::abs(Vofx[lLowI - 1 - i]);
            Ihigh[i] = std::abs(Iofx[lLowI - 1 - i]);
        }
        for (size_t i = 0; i < lLengthLow; i++) {
            Vlow[i] = std::abs(Vofx[lLowI + i]);
            Ilow[i] = std::abs(Iofx[lLowI + i]);
        }
    }
    std::valarray<double> Irex(lLengthLow), Vrex(lLengthLow);
    writeIV("test1.txt", Ilow, Vlow);
    writeIV("test2.txt", Ihigh, Vhigh);
    Resample(Vrex, Irex, Vhigh, Ihigh, Vlow, Ilow);
    writeIV("test3.txt", Irex, Vrex);
    return ChiSqHi(Vlow, Vrex);
}
