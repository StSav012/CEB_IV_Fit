#ifndef MINS_H
#define MINS_H

#include <functional>

class Bracketmethod {
public:
    double ax, bx, cx;

    void setABC(double a, double b, double c);

    void bracket(double a, double b, const std::function<double(double)>& f);
};

std::tuple<double, double> GoldenMinimize(const std::function<double(double)>& f,
                                          double limit1, double limit2,
                                          double initialGuess,
                                          double tolerance = 3.0e-8);

std::tuple<double, double> BrentMinimize(const std::function<double(double)>& f,
                                         double limit1, double limit2,
                                         double initialGuess,
                                         double tolerance = 3.0e-8);

std::tuple<double, double> DBrentMinimize(const std::function<double(double)>& f,
                                          const std::function<double(double)>& df,
                                          double limit1, double limit2,
                                          double initialGuess,
                                          double tolerance = 3.0e-8);

#endif
