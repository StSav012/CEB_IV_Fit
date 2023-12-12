#ifndef CFOO_H
#define CFOO_H

#include <valarray>

class CFoo {
    const size_t parCount = 21;

    std::valarray<double> Iexp, Vexp;
    std::valarray<double> Inum, Vnum;

    std::valarray<double> par;
    std::valarray<bool> ToFit;
    size_t iParNum;

public:
    double operator()(double dParam);

    size_t CEB_2eq_parallel_lite();

    void SeqFit(size_t runCount, const std::valarray<double>& Irex);

    size_t loadExperimentData(const std::string& filename, bool removeOffset = false);

    [[nodiscard]] std::tuple<std::valarray<double>, std::valarray<double>> resample() const;

    explicit CFoo(size_t parIndex);
};

#endif
