#include <cmath>

#include "nrutil_nr.h"

struct Bracketmethod
{
    double ax, bx, cx, fa, fb, fc;
    template<class T>
    void bracket(const double a, const double b, T &func)
    {
        const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
        ax = a;
        bx = b;
        double fu;
        fa = func(ax);
        fb = func(bx);
        if (fb > fa) {
            std::swap(ax, bx);
            std::swap(fb, fa);
        }
        cx = bx + GOLD * (bx - ax);
        fc = func(cx);
        while (fb > fc) {
            double r = (bx - ax) * (fb - fc);
            double q = (bx - cx) * (fb - fa);
            double u = bx
                       - ((bx - cx) * q - (bx - ax) * r)
                             / (2.0 * SIGN(std::max(std::abs(q - r), TINY), q - r));
            double ulim = bx + GLIMIT * (cx - bx);
            if ((bx - u) * (u - cx) > 0.0) {
                fu = func(u);
                if (fu < fc) {
                    ax = bx;
                    bx = u;
                    fa = fb;
                    fb = fu;
                    return;
                } else if (fu > fb) {
                    cx = u;
                    fc = fu;
                    return;
                }
                u = cx + GOLD * (cx - bx);
                fu = func(u);
            } else if ((cx - u) * (u - ulim) > 0.0) {
                fu = func(u);
                if (fu < fc) {
                    shft3(bx, cx, u, u + GOLD * (u - cx));
                    shft3(fb, fc, fu, func(u));
                }
            } else if ((u - ulim) * (ulim - cx) >= 0.0) {
                u = ulim;
                fu = func(u);
            } else {
                u = cx + GOLD * (cx - bx);
                fu = func(u);
            }
            shft3(ax, bx, cx, u);
            shft3(fa, fb, fc, fu);
        }
    }
    inline void shft2(double &a, double &b, const double c)
    {
        a = b;
        b = c;
    }
    inline void shft3(double &a, double &b, double &c, const double d)
    {
        a = b;
        b = c;
        c = d;
    }
    inline void mov3(double &a, double &b, double &c, const double d, const double e, const double f)
    {
        a = d;
        b = e;
        c = f;
    }
};

struct Golden : Bracketmethod
{
    double xmin, fmin;
    const double tol;
    Golden(const double toll = 3.0e-8)
        : tol(toll)
    {}
    template<class T>
    double minimize(T &func)
    {
        clock_t start, finish;
        start = clock();
        const double R = 0.61803399, C = 1.0 - R;
        double x1, x2;
        double x0 = ax;
        double x3 = cx;
        if (std::abs(cx - bx) > std::abs(bx - ax)) {
            x1 = bx;
            x2 = bx + C * (cx - bx);
        } else {
            x2 = bx;
            x1 = bx - C * (bx - ax);
        }
        double f1 = func(x1);
        double f2 = func(x2);
        while (std::abs(x3 - x0) > tol * (std::abs(x1) + std::abs(x2))) {
            if (f2 < f1) {
                shft3(x0, x1, x2, R * x2 + C * x3);
                shft2(f1, f2, func(x2));
                finish = clock();
                printf("\nCURRENT XMIN:%lf CHISQMIN: %e\n", x2, f2);
                FILE *conv = fopen("converg.txt", "a");
                fprintf(conv, "%lf %e %lf\n", x2, f2, (double) (finish - start) / CLOCKS_PER_SEC);
                fclose(conv);
                //getch();
            } else {
                shft3(x3, x2, x1, R * x1 + C * x0);
                shft2(f2, f1, func(x1));
                finish = clock();
                printf("\nCURRENT XMIN:%lf CHISQMIN: %e\n", x1, f1);
                FILE *conv = fopen("converg.txt", "a");
                fprintf(conv, "%lf %e %lf\n", x1, f1, (double) (finish - start) / CLOCKS_PER_SEC);
                fclose(conv);
                //getch();
            }
        }
        if (f1 < f2) {
            xmin = x1;
            fmin = f1;
        } else {
            xmin = x2;
            fmin = f2;
        }
        return xmin;
    }
};

struct Brent : Bracketmethod
{
    double xmin, fmin;
    const double tol;
    Brent(const double toll = 3.0e-8)
        : tol(toll)
    {}
    template<class T>
    double minimize(T &func)
    {
        const int ITMAX = 100;
        const double CGOLD = 0.3819660;
        const double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
        double a, b, d = 0.0, etemp, fu, fv, fw, fx;
        double p, q, r, tol1, tol2, u, v, w, x, xm;
        double e = 0.0;

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);
        x = w = v = bx;
        fw = fv = fx = func(x);
        clock_t start, finish;
        start = clock();
        for (int iter = 0; iter < ITMAX; iter++) {
            xm = 0.5 * (a + b);
            tol2 = 2.0 * (tol1 = tol * std::abs(x) + ZEPS);
            if (std::abs((x - xm) / xm) <= tol) {
                fmin = fx;
                return xmin = x;
            }
            if (std::abs(e) > tol1) {
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0)
                    p = -p;
                q = std::abs(q);
                etemp = e;
                e = d;
                if (std::abs(p) >= std::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                    d = CGOLD * (e = (x >= xm ? a - x : b - x));
                else {
                    d = p / q;
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2)
                        d = SIGN(tol1, xm - x);
                }
            } else {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            }
            u = (std::abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
            fu = func(u);
            if (fu <= fx) {
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
            finish = clock();
            printf("\nCURRENT XMIN:%lf CHISQMIN: %e\n", x, fx);
            FILE *conv = fopen("converg.txt", "a");
            fprintf(conv, "%lf %e %lf\n", x, fx, (double) (finish - start) / CLOCKS_PER_SEC);
            fclose(conv);
        }
        throw("Too many iterations in brent");
    }
};

struct Dbrent : Bracketmethod
{
    double xmin, fmin;
    const double tol;
    Dbrent(const double toll = 3.0e-8)
        : tol(toll)
    {}
    template<class T>
    double minimize(T &funcd)
    {
        const int ITMAX = 100;
        const double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
        bool ok1, ok2;
        double a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
        double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);
        x = w = v = bx;
        fw = fv = fx = funcd(x);
        dw = dv = dx = funcd.df(x);
        for (int iter = 0; iter < ITMAX; iter++) {
            xm = 0.5 * (a + b);
            tol1 = tol * std::abs(x) + ZEPS;
            tol2 = 2.0 * tol1;
            if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                fmin = fx;
                return xmin = x;
            }
            if (std::abs(e) > tol1) {
                d1 = 2.0 * (b - a);
                d2 = d1;
                if (dw != dx)
                    d1 = (w - x) * dx / (dx - dw);
                if (dv != dx)
                    d2 = (v - x) * dx / (dx - dv);
                u1 = x + d1;
                u2 = x + d2;
                ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
                ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
                olde = e;
                e = d;
                if (ok1 || ok2) {
                    if (ok1 && ok2)
                        d = (std::abs(d1) < std::abs(d2) ? d1 : d2);
                    else if (ok1)
                        d = d1;
                    else
                        d = d2;
                    if (std::abs(d) <= std::abs(0.5 * olde)) {
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                            d = SIGN(tol1, xm - x);
                    } else {
                        d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                    }
                } else {
                    d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                }
            } else {
                d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
            }
            if (std::abs(d) >= tol1) {
                u = x + d;
                fu = funcd(u);
            } else {
                u = x + SIGN(tol1, d);
                fu = funcd(u);
                if (fu > fx) {
                    fmin = fx;
                    return xmin = x;
                }
            }
            du = funcd.df(u);
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
        throw("Too many iterations in routine dbrent");
    }
};
