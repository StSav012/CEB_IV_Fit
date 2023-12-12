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

double COff::operator()(const double offset) const {
    /*
     *  Make an object of the class act like a callable
     *  
     *  TODO: describe what the call is doing.
     */
    const size_t count = Iref.size();
    if ((Iref[std::slice(0, count - 1, 1)] >= Iref[std::slice(1, count - 1, 1)]).sum()) {
        throw std::runtime_error("I doesn't rise monotonously");
    }

    const std::valarray<double> Vofx = Vref - offset;
    const std::valarray<double> Iofx = Iref;

    // find the index of the item closest to zero, but still positive
    const std::valarray<double> absIofx = std::abs(Iofx);
    size_t indexOfTheSmallestAbsIofx = 0;
    double smallestAbsIofx = absIofx[indexOfTheSmallestAbsIofx];
    for (auto pAbsIofxItem = begin(absIofx); pAbsIofxItem != end(absIofx); ++pAbsIofxItem) {
        if (*pAbsIofxItem >= 0 && *pAbsIofxItem < smallestAbsIofx) {
            indexOfTheSmallestAbsIofx = std::distance(pAbsIofxItem, begin(absIofx));
            smallestAbsIofx = *pAbsIofxItem;
        }
    }

    // split the IV curve into positive and negative branches
    const size_t lLengthPos = count - indexOfTheSmallestAbsIofx;
    const size_t lLengthNeg = indexOfTheSmallestAbsIofx;
    std::valarray<double> Vlow, Ilow, Vhigh, Ihigh;
    size_t lLengthLow, lLengthHigh;
    if (lLengthPos > lLengthNeg) {
        // i.e., indexOfTheSmallestAbsIofx < count / 2
        lLengthHigh = lLengthPos;
        lLengthLow = lLengthNeg;
        Vhigh = Vofx[std::slice(indexOfTheSmallestAbsIofx, lLengthHigh, 1)];
        Ihigh = Iofx[std::slice(indexOfTheSmallestAbsIofx, lLengthHigh, 1)];
        Vlow = Vofx[std::slice(indexOfTheSmallestAbsIofx - 1, lLengthLow, -1)];
        Ilow = Iofx[std::slice(indexOfTheSmallestAbsIofx - 1, lLengthLow, -1)];
    } else {
        lLengthHigh = lLengthNeg;
        lLengthLow = lLengthPos;
        Vhigh = Vofx[std::slice(indexOfTheSmallestAbsIofx - 1, lLengthHigh, -1)];
        Ihigh = Iofx[std::slice(indexOfTheSmallestAbsIofx - 1, lLengthHigh, -1)];
        Vlow = Vofx[std::slice(indexOfTheSmallestAbsIofx, lLengthLow, 1)];
        Ilow = Iofx[std::slice(indexOfTheSmallestAbsIofx, lLengthLow, 1)];
    }

    writeIV("IV low.txt", Ilow, Vlow);
    writeIV("IV high.txt", Ihigh, Vhigh);

    // symmetrize the voltage points of the branches
    auto [Irex, Vrex] = Resample(Ihigh, Vhigh, Ilow, Vlow);
    writeIV("IV rex.txt", Irex, Vrex);

    // compare the shortest branch with the resampled one
    return ChiSqHi(Ilow, Irex);
}
