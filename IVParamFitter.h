#ifndef CFOO_H
#define CFOO_H

#include <string>
#include <unordered_map>
#include <valarray>

class IVParamFitter {
    std::valarray<double> Iexp, Vexp;
    std::valarray<double> Inum, Vnum;

    std::unordered_map<std::string, double> par;
    std::unordered_map<std::string, bool> ToFit;
    std::string parameterName;

public:
    double operator()(double dParam);

    size_t computeCEBProperties();

    void SeqFit(size_t runCount, const std::valarray<double>& Irex);

    size_t loadExperimentData(const std::string& filename, bool removeOffset = false);

    [[nodiscard]] std::tuple<std::valarray<double>, std::valarray<double>> resample() const;

    explicit IVParamFitter();
};

#endif
