#pragma once
#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include "constants.hpp"
#include <cmath>
#include <iostream>

// some fast math algorathms used in the code.
namespace basic_math
{
constexpr unsigned long long factorial_table[] = {1,
                                                  1,
                                                  2,
                                                  6,
                                                  24,
                                                  120,
                                                  720,
                                                  5040,
                                                  40320,
                                                  362880,
                                                  3628800,
                                                  39916800,
                                                  479001600,
                                                  6227020800,
                                                  87178291200,
                                                  1307674368000,
                                                  20922789888000,
                                                  355687428096000,
                                                  6402373705728000,
                                                  121645100408832000,
                                                  2432902008176640000};

constexpr unsigned long long doublefactorial_table[] = {1,
                                                        1,
                                                        2,
                                                        3,
                                                        8,
                                                        15,
                                                        48,
                                                        105,
                                                        384,
                                                        945,
                                                        3840,
                                                        10395,
                                                        46080,
                                                        135135,
                                                        645120,
                                                        2027025,
                                                        10321920,
                                                        34459425,
                                                        185794560,
                                                        654729075,
                                                        3715891200,
                                                        13749310575,
                                                        81749606400,
                                                        316234143225,
                                                        1961990553600,
                                                        7905853580625,
                                                        51011754393600,
                                                        213458046676875,
                                                        1428329123020800,
                                                        6190283353629375,
                                                        42849873690624000};

constexpr double precomputed_sqrt[] = {1.0,
                                       1.41421356237309504880,
                                       1.73205080756887729352,
                                       2.0,
                                       2.23606797749978969640,
                                       2.44948974278317809819,
                                       2.64575131106459059050,
                                       2.82842712474619009760,
                                       3.0,
                                       3.16227766016837933199,
                                       3.31662479035539984911,
                                       3.46410161513775458705,
                                       3.60555127546398929311,
                                       3.74165738677394138558,
                                       3.87298334620741688517,
                                       4.0,
                                       4.12310562561766054982,
                                       4.24264068711928514640,
                                       4.35889894354067355223,
                                       4.47213595499957939281,
                                       4.58257569495584000658,
                                       4.69041575982342955456};

//* quick factorial function n!, quick for n<=20.
unsigned long long factorial(unsigned int n)
{
    if (n <= 20)
    {
        return factorial_table[n];
    }
    return n * factorial(n - 1);
}

//* quick double-factorial function n!!, quick for n<=30.
unsigned long long doublefactorial(unsigned int n)
{
    if (n <= 30)
    {
        return doublefactorial_table[n];
    }
    return n * doublefactorial(n - 2);
}

//* quick pow(2,n) function
unsigned long long pow_two(unsigned int n) { return (1ULL << n); }

//* quick pow(x,n) function
double quick_pow(double x, int n)
{
    bool neg = (n < 0);
    double ans = 1.0;
    n = std::abs(n);
    while (n)
    {
        if (n & 1)
        {
            ans *= x;
        }
        n >>= 1;
        x *= x;
    }
    return neg ? 1.0 / ans : ans;
}

// fast pow function pow(a,b), needed to be checked.
inline double fastPow(double a, double b)
{
    union
    {
        double d;
        int x[2];
    } u = {a};
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
}

// integer from double.
// --------------------
int make_int(const double x_int) { return static_cast<int>(rint(x_int)); }

// -1 pow integer
// --------------
int minus_one_pow(const int x_int) { return ((x_int & 0x01) ? (-1) : (1)); }

double minus_one_pow_double(const int x_int) { return ((x_int & 0x01) ? (-1.0) : (1.0)); }

// return (-1)^n
double iphase(const int n) { return ((n & 0x01) ? (-1.0) : (1.0)); }

// sign function
// -------------
int SIGN(const double x) { return ((x < 0) ? (-1) : (1)); }

// unsigned integer from double.
// ---------------------------- -
unsigned int make_uns_int(const double x_int) { return static_cast<unsigned int>(rint(x_int)); }

inline double delta_double(double m, double n) { return (std::abs(m - n) < 1E-5) ? 1.0 : 0.0; }

//* Quick calculation of hat(J) = sqrt[2.J + 1],
//* It is quick for J <= 21/2.
double hat(const double j)
{
    const int two_j = static_cast<int>(2 * j);
    if (two_j <= 21)
    {
        return precomputed_sqrt[two_j];
    }
    else
    {
        return sqrt(2 * j + 1);
    }
}

// log Gamma function.
//---------------------------------------------------------------
std::complex<double> log_Gamma(const std::complex<double> &z)
{
    const double x = real(z);
    const double y = imag(z);

    if (x >= 0.5)
    {
        const double log_sqrt_2Pi = 0.91893853320467274177;

        const double g = 4.7421875;

        const std::complex<double> z_m_0p5 = z - 0.5;

        const std::complex<double> z_pg_m0p5 = z_m_0p5 + g;

        const std::complex<double> zm1 = z - 1.0;

        const double c[15] = {0.99999999999999709182,     57.156235665862923517,      -59.597960355475491248,
                              14.136097974741747174,      -0.49191381609762019978,    0.33994649984811888699e-4,
                              0.46523628927048575665e-4,  -0.98374475304879564677e-4, 0.15808870322491248884e-3,
                              -0.21026444172410488319e-3, 0.21743961811521264320e-3,  -0.16431810653676389022e-3,
                              0.84418223983852743293e-4,  -0.26190838401581408670e-4, 0.36899182659531622704e-5};

        std::complex<double> sum = c[0];

        for (int i = 1; i < 15; i++)
            sum += c[i] / (zm1 + double(i));

        const std::complex<double> log_Gamma_z = log_sqrt_2Pi + log(sum) + z_m_0p5 * log(z_pg_m0p5) - z_pg_m0p5;

        return log_Gamma_z;
    }
    else if (y >= 0.0)
    {
        const int n = (x < rint(x)) ? (static_cast<int>(rint(x)) - 1) : (static_cast<int>(rint(x)));

        const double log_Pi = 1.1447298858494002;

        const std::complex<double> log_const(-0.693147180559945309417, constants::pi / 2.0);

        const std::complex<double> i_Pi(0.0, constants::pi);

        const std::complex<double> eps = z - double(n);

        const std::complex<double> log_sin_Pi_z =
            (y > 110) ? (-i_Pi * z + log_const) : (log(sin(constants::pi * eps)) - i_Pi * double(n));

        const std::complex<double> log_Gamma_z = log_Pi - log_sin_Pi_z - basic_math::log_Gamma(1.0 - z);

        return log_Gamma_z;
    }
    else
    {
        return conj(basic_math::log_Gamma(conj(z)));
    }
}
//---------------------------------------------------------------

// Gamma function.
std::complex<double> Gamma(const std::complex<double> &z) { return exp(basic_math::log_Gamma(z)); }

// Coulomb phase shift.
// --------------------
// It is given by the formula [Gamma[1 + l + I.eta] - Gamma[1 + l - I.eta]]/[2i].
// 0 is returned if 1 + l +/- I.eta is a negative integer.
//
// Variables:
// ----------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta , one_over_two_I : i.eta , 1/[2.i] .
// arg_plus , arg_minus : 1 + l + i.eta , 1 + l - i.eta.
// log_Gamma_plus , log_Gamma_minus : logs of Gamma[1 + l + I.eta] , Gamma[1 + l - I.eta].
// sigma_l : returned result.
//---------------------------------------------------------------
std::complex<double> sigma_l_calc(const int &l, double nc)
{
    std::complex<double> l_complex(l, 0);
    std::complex<double> eta(nc, 0);
    const std::complex<double> Ieta(-imag(eta), real(eta));
    const std::complex<double> one_over_two_I(0, -0.5);

    const std::complex<double> arg_plus = 1.0 + l + Ieta;
    const std::complex<double> arg_minus = 1.0 + l - Ieta;

    if ((rint(real(arg_plus)) == arg_plus) && (rint(real(arg_plus)) <= 0.0))
    {
        return 0.0;
    }

    if ((rint(real(arg_minus)) == arg_minus) && (rint(real(arg_minus)) <= 0.0))
    {
        return 0.0;
    }

    const std::complex<double> log_Gamma_plus = basic_math::log_Gamma(arg_plus);
    const std::complex<double> log_Gamma_minus = basic_math::log_Gamma(arg_minus);

    const std::complex<double> sigma_l = (log_Gamma_plus - log_Gamma_minus) * one_over_two_I;

    return sigma_l;
}
//---------------------------------------------------------------

} // namespace basic_math

#endif // BASIC_MATH_HPP