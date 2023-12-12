#include <cmath>
#include <chrono>
#include <fstream>
#include <functional>

#include "mins.h"
#include "utils.h"

static void shft2(double& a, double& b, const double c) {
    a = b;
    b = c;
}

static void shft3(double& a, double& b, double& c, const double d) {
    a = b;
    b = c;
    c = d;
}

static void mov3(double& a, double& b, double& c, const double new_a, const double new_b, const double new_c) {
    a = new_a;
    b = new_b;
    c = new_c;
}

void Bracketmethod::setABC(const double a, const double b, const double c) {
    ax = a;
    bx = b;
    cx = c;
}

void Bracketmethod::bracket(const double a, const double b, const std::function<double(double)>& f) {
    const double GOLD = (std::sqrt(5.0) + 1.0) / 2.0;
    ax = a;
    bx = b;
    double fu;
    double fa = f(ax);
    double fb = f(bx);
    if (fb > fa) {
        std::swap(ax, bx);
        std::swap(fb, fa);
    }
    cx = bx + GOLD * (bx - ax);
    double fc = f(cx);
    while (fb > fc) {
        constexpr double TINY = 1.0e-20;
        constexpr double GLIMIT = 100.0;
        const double r = (bx - ax) * (fb - fc);
        const double q = (bx - cx) * (fb - fa);
        double u = bx
                   - ((bx - cx) * q - (bx - ax) * r)
                   / (2.0 * std::copysign(std::max(std::abs(q - r), TINY), q - r));
        double ulim = bx + GLIMIT * (cx - bx);
        if ((bx > u) == (u > cx)) {
            fu = f(u);
            if (fu < fc) {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            }
            if (fu > fb) {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + GOLD * (cx - bx);
            fu = f(u);
        } else if ((cx > u) == (u > ulim)) {
            fu = f(u);
            if (fu < fc) {
                shft3(bx, cx, u, u + GOLD * (u - cx));
                shft3(fb, fc, fu, f(u));
            }
        } else if ((u >= ulim) == (ulim >= cx)) {
            u = ulim;
            fu = f(u);
        } else {
            u = cx + GOLD * (cx - bx);
            fu = f(u);
        }
        shft3(ax, bx, cx, u);
        shft3(fa, fb, fc, fu);
    }
}

std::tuple<double, double> GoldenMinimize(const std::function<double(double)>& f,
                                          const double limit1, const double limit2,
                                          const double initialGuess, const double tolerance) {
    const std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    const double R = (std::sqrt(5.0) - 1.0) / 2.0, C = 1.0 - R;
    double x1, x2;
    double x0 = limit1;
    double x3 = limit2;
    if (std::abs(limit2 - initialGuess) > std::abs(initialGuess - limit1)) {
        x1 = initialGuess;
        x2 = initialGuess + C * (limit2 - initialGuess);
    } else {
        x2 = initialGuess;
        x1 = initialGuess - C * (initialGuess - limit1);
    }
    double f1 = f(x1);
    double f2 = f(x2);
    while (std::abs(x3 - x0) > tolerance * (std::abs(x1) + std::abs(x2))) {
        if (f2 < f1) {
            shft3(x0, x1, x2, R * x2 + C * x3);
            shft2(f1, f2, f(x2));
            writeConverg(x2, f2, start);
            // std::getchar();
        } else {
            shft3(x3, x2, x1, R * x1 + C * x0);
            shft2(f2, f1, f(x1));
            writeConverg(x1, f1, start);
            // std::getchar();
        }
    }
    if (f1 < f2) {
        return std::make_tuple(x1, f1);
    }
    return std::make_tuple(x2, f2);
}

std::tuple<double, double> BrentMinimize(const std::function<double(double)>& f,
                                         const double limit1, const double limit2,
                                         const double initialGuess, const double tolerance) {
    constexpr int ITMAX = 100;
    constexpr double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
    const double C = 1.0 - (std::sqrt(5.0) - 1.0) / 2.0;

    double d = 0.0, fv, fx;
    double tol1, u, v, w;
    double e = 0.0;

    double a = std::min(limit1, limit2);
    double b = std::max(limit1, limit2);
    double x = w = v = initialGuess;
    double fw = fv = fx = f(x);
    const std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    for (int iter = 0; iter < ITMAX; ++iter) {
        const double xm = 0.5 * (a + b);
        double tol2 = 2.0 * (tol1 = tolerance * std::abs(x) + ZEPS);
        if (std::abs((x - xm) / xm) <= tolerance) {
            return std::make_tuple(x, fx);
        }
        if (std::abs(e) > tol1) {
            const double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = std::abs(q);
            const double etemp = e;
            e = d;
            if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                e = x >= xm ? a - x : b - x;
                d = C * e;
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = std::copysign(tol1, xm - x);
            }
        } else {
            e = x >= xm ? a - x : b - x;
            d = C * e;
        }
        u = (std::abs(d) >= tol1 ? x + d : x + std::copysign(tol1, d));
        if (const double fu = f(u); fu <= fx) {
            if (u >= x)
                a = x;
            else
                b = x;
            shft3(v, w, x, u);
            shft3(fv, fw, fx, fu);
        } else {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
        writeConverg(x, fx, start);
    }
    throw std::runtime_error("Too many iterations in brent");
}

std::tuple<double, double> DBrentMinimize(const std::function<double(double)>& f,
                                          const std::function<double(double)>& df,
                                          const double limit1, const double limit2,
                                          const double initialGuess, const double tolerance) {
    constexpr int ITMAX = 100;
    constexpr double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
    double d = 0.0, dv, dx, e = 0.0;
    double fu, fv, fx, u, v, w;

    double a = std::min(limit1, limit2);
    double b = std::max(limit1, limit2);
    double x = w = v = initialGuess;
    double fw = fv = fx = f(x);
    double dw = dv = dx = df(x);
    for (int iter = 0; iter < ITMAX; iter++) {
        const double xm = 0.5 * (a + b);
        double tol1 = tolerance * std::abs(x) + ZEPS;
        double tol2 = 2.0 * tol1;
        if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            return std::make_tuple(x, fx);
        }
        if (std::abs(e) > tol1) {
            double d1 = 2.0 * (b - a);
            double d2 = d1;
            if (dw != dx)
                d1 = (w - x) * dx / (dx - dw);
            if (dv != dx)
                d2 = (v - x) * dx / (dx - dv);
            const double u1 = x + d1;
            const double u2 = x + d2;
            const bool ok1 = (a > u1) == (u1 > b) && std::signbit(dx) != std::signbit(d1);
            const bool ok2 = (a > u2) == (u2 > b) && std::signbit(dx) != std::signbit(d2);
            const double olde = e;
            e = d;
            if (ok1 || ok2) {
                if (ok1 && ok2)
                    d = std::abs(d1) < std::abs(d2) ? d1 : d2;
                else if (ok1)
                    d = d1;
                else
                    d = d2;
                if (std::abs(d) <= std::abs(0.5 * olde)) {
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2)
                        d = std::copysign(tol1, xm - x);
                } else {
                    e = dx >= 0.0 ? a - x : b - x;
                    d = 0.5 * e;
                }
            } else {
                e = dx >= 0.0 ? a - x : b - x;
                d = 0.5 * e;
            }
        } else {
            e = dx >= 0.0 ? a - x : b - x;
            d = 0.5 * e;
        }
        if (std::abs(d) >= tol1) {
            u = x + d;
            fu = f(u);
        } else {
            u = x + std::copysign(tol1, d);
            fu = f(u);
            if (fu > fx) {
                return std::make_tuple(x, fx);
            }
        }
        const double du = df(u);
        if (fu <= fx) {
            if (u >= x)
                a = x;
            else
                b = x;
            mov3(v, fv, dv, w, fw, dw);
            mov3(w, fw, dw, x, fx, dx);
            mov3(x, fx, dx, u, fu, du);
        } else {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x) {
                mov3(v, fv, dv, w, fw, dw);
                mov3(w, fw, dw, u, fu, du);
            } else if (fu < fv || v == x || v == w) {
                mov3(v, fv, dv, u, fu, du);
            }
        }
    }
    throw std::runtime_error("Too many iterations in routine dbrent");
}
