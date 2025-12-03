#pragma once
#ifndef INTERACTION_aPWD_HPP
#define INTERACTION_aPWD_HPP

#include "lib_define.hpp"
#include "mymath.hpp"
#include "state.hpp"

namespace interaction_aPWD
{
using std::string;

// Error message written and all processes stop
void error_message_print_abort(const string &error_message) { std::cout << error_message << std::endl; }

constexpr double ppi = 3.14159265358979323846;

// partial-wave projection by aPWD method under SYM-LSJ basis convention.
double V_kernal_SYM(const states::LSJ_State &s1, const states::LSJ_State &s2, const double &f1, const double &f2,
                    const double &f3, const double &f4, const double &f5, const double &f6, const double &f7,
                    const double &f8, const double &x)
{
    const double pmag = s1.momentum;
    const double ppmag = s2.momentum;

    if (!(s1.j == s2.j && s1.s == s2.s && s1.j >= 0 && s2.j >= 0 && s1.l >= 0 && s2.l >= 0 && s1.s >= 0 && s2.s >= 0 &&
          std::abs(s1.l - s2.l) <= 2))
    {
        return 0;
    }
    if (s2.l == 0 && s1.l == 0 && s1.s == 0 && s1.j == 0)
    {
        return -2 * ppi *
               (-f1 + 3 * f2 + pow(pmag, 2) * pow(ppmag, 2) * f4 - pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 +
                pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 + 2 * pmag * ppmag * x * f5 + pow(pmag, 2) * f6 +
                pow(ppmag, 2) * f6 - 2 * pmag * ppmag * x * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7);
    }
    if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 0)
    {
        return -2 * ppi *
               (-2 * pmag * ppmag * pow(x, 2) * f3 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 3) * f4 +
                2 * pmag * ppmag * (f3 + f5 - f6) +
                x * (-f1 - f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 +
                     pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7));
    }
    if (s2.l == 1 && s1.l == 1 && s1.s == 0 && s1.j == 1)
    {
        return 2 * ppi * x *
               (f1 - 3 * f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                pow(pmag, 2) * f5 - pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - pow(pmag, 2) * f6 -
                pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6 + 2 * pow(pmag, 2) * f7 - 2 * pow(ppmag, 2) * f7);
    }
    if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 1)
    {
        return 2 * ppi *
               (pmag * ppmag * pow(x, 2) * (f3 + f5 - f6) - pmag * ppmag * (f3 - f5 + f6) +
                x * (f1 + f2 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 + pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 -
                     2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7));
    }
    if (s2.l == 0 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return (2 * ppi *
                (3 * f1 + 3 * f2 + pow(pmag, 2) * pow(ppmag, 2) * f4 - pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 +
                 pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 + 2 * pmag * ppmag * x * f5 + pow(pmag, 2) * f6 +
                 pow(ppmag, 2) * f6 - 2 * pmag * ppmag * x * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7)) /
               3.;
    }
    if (s2.l == 0 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return -(2 * sqrt(2) * ppi *
                 (4 * pmag * ppmag * x * (f5 - f6) +
                  pow(pmag, 2) * (pow(ppmag, 2) * (-1 + pow(x, 2)) * f4 + 2 * (f5 + f6 - 2 * f7)) +
                  pow(ppmag, 2) * (-1 + 3 * pow(x, 2)) * (f5 + f6 + 2 * f7))) /
               3.;
    }
    if (s2.l == 2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
    {
        return -(2 * sqrt(2) * ppi *
                 (4 * pmag * ppmag * x * (f5 - f6) +
                  pow(pmag, 2) * (pow(ppmag, 2) * (-1 + pow(x, 2)) * f4 + (-1 + 3 * pow(x, 2)) * (f5 + f6 - 2 * f7)) +
                  2 * pow(ppmag, 2) * (f5 + f6 + 2 * f7))) /
               3.;
    }
    if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 1)
    {
        return (ppi *
                ((-3 + 9 * pow(x, 2)) * f1 + (-3 + 9 * pow(x, 2)) * f2 - 18 * pmag * ppmag * x * f3 +
                 18 * pmag * ppmag * pow(x, 3) * f3 - 5 * pow(pmag, 2) * pow(ppmag, 2) * f4 +
                 14 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                 9 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 4) * f4 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 -
                 4 * pmag * ppmag * x * f5 - 3 * pow(pmag, 2) * pow(x, 2) * f5 - 3 * pow(ppmag, 2) * pow(x, 2) * f5 +
                 pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 + 4 * pmag * ppmag * x * f6 -
                 3 * pow(pmag, 2) * pow(x, 2) * f6 - 3 * pow(ppmag, 2) * pow(x, 2) * f6 - 2 * pow(pmag, 2) * f7 +
                 2 * pow(ppmag, 2) * f7 + 6 * pow(pmag, 2) * pow(x, 2) * f7 - 6 * pow(ppmag, 2) * pow(x, 2) * f7)) /
               3.;
    }
    if (s2.l == 2 && s1.l == 2 && s1.s == 0 && s1.j == 2)
    {
        return ppi * (-1 + 3 * pow(x, 2)) *
               (f1 - 3 * f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                pow(pmag, 2) * f5 - pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - pow(pmag, 2) * f6 -
                pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6 + 2 * pow(pmag, 2) * f7 - 2 * pow(ppmag, 2) * f7);
    }
    if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 2)
    {
        return ppi *
               ((-1 + 3 * pow(x, 2)) * f1 + (-1 + 3 * pow(x, 2)) * f2 - 2 * pmag * ppmag * x * f3 +
                2 * pmag * ppmag * pow(x, 3) * f3 + pow(pmag, 2) * pow(ppmag, 2) * f4 -
                2 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 4) * f4 -
                pow(pmag, 2) * f5 - pow(ppmag, 2) * f5 + 3 * pow(pmag, 2) * pow(x, 2) * f5 +
                3 * pow(ppmag, 2) * pow(x, 2) * f5 + 4 * pmag * ppmag * pow(x, 3) * f5 - pow(pmag, 2) * f6 -
                pow(ppmag, 2) * f6 + 3 * pow(pmag, 2) * pow(x, 2) * f6 + 3 * pow(ppmag, 2) * pow(x, 2) * f6 -
                4 * pmag * ppmag * pow(x, 3) * f6 + 2 * pow(pmag, 2) * f7 - 2 * pow(ppmag, 2) * f7 -
                6 * pow(pmag, 2) * pow(x, 2) * f7 + 6 * pow(ppmag, 2) * pow(x, 2) * f7);
    }
    if (s2.l == 1 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return (-2 * ppi *
                (2 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 3) * f4 + pmag * ppmag * (-5 * f3 + f5 - f6) +
                 pmag * ppmag * pow(x, 2) * (5 * f3 - 3 * f5 + 3 * f6) -
                 x * (5 * f1 + 5 * f2 + 2 * pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 +
                      pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7))) /
               5.;
    }
    if (s2.l == 1 && s1.l == 3 && s1.s == 1 && s1.j == 2)
    {
        return -(2 * sqrt(6) * ppi *
                 (2 * pmag * ppmag * (-1 + 3 * pow(x, 2)) * (f5 - f6) +
                  pow(pmag, 2) * x * (pow(ppmag, 2) * (-1 + pow(x, 2)) * f4 + 2 * (f5 + f6 - 2 * f7)) +
                  pow(ppmag, 2) * x * (-3 + 5 * pow(x, 2)) * (f5 + f6 + 2 * f7))) /
               5.;
    }
    if (s2.l == 3 && s1.l == 1 && s1.s == 1 && s1.j == 2)
    {
        return -(2 * sqrt(6) * ppi *
                 (2 * pmag * ppmag * (-1 + 3 * pow(x, 2)) * (f5 - f6) +
                  pow(pmag, 2) * x *
                      (pow(ppmag, 2) * (-1 + pow(x, 2)) * f4 + (-3 + 5 * pow(x, 2)) * (f5 + f6 - 2 * f7)) +
                  2 * pow(ppmag, 2) * x * (f5 + f6 + 2 * f7))) /
               5.;
    }
    if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 2)
    {
        return -0.2 *
               (ppi * (-50 * pmag * ppmag * pow(x, 4) * f3 + 25 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 5) * f4 -
                       2 * pmag * ppmag * (5 * f3 + f5 - f6) + 6 * pmag * ppmag * pow(x, 2) * (10 * f3 + f5 - f6) +
                       x * (15 * f1 + 15 * f2 + 19 * pow(pmag, 2) * pow(ppmag, 2) * f4 - 3 * pow(pmag, 2) * f5 -
                            3 * pow(ppmag, 2) * f5 - 3 * pow(pmag, 2) * f6 - 3 * pow(ppmag, 2) * f6 +
                            6 * pow(pmag, 2) * f7 - 6 * pow(ppmag, 2) * f7) +
                       pow(x, 3) * (-25 * f1 - 25 * f2 - 44 * pow(pmag, 2) * pow(ppmag, 2) * f4 +
                                    5 * pow(pmag, 2) * f5 + 5 * pow(ppmag, 2) * f5 + 5 * pow(pmag, 2) * f6 +
                                    5 * pow(ppmag, 2) * f6 - 10 * pow(pmag, 2) * f7 + 10 * pow(ppmag, 2) * f7)));
    }
    if (s2.l == 3 && s1.l == 3 && s1.s == 0 && s1.j == 3)
    {
        return ppi * x * (-3 + 5 * pow(x, 2)) *
               (f1 - 3 * f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                pow(pmag, 2) * f5 - pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - pow(pmag, 2) * f6 -
                pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6 + 2 * pow(pmag, 2) * f7 - 2 * pow(ppmag, 2) * f7);
    }
    if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 3)
    {
        return (ppi * (5 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 5) * f4 +
                       5 * pmag * ppmag * pow(x, 4) * (f3 + 3 * f5 - 3 * f6) -
                       6 * pmag * ppmag * pow(x, 2) * (f3 + f5 - f6) + pmag * ppmag * (f3 - f5 + f6) +
                       10 * pow(x, 3) *
                           (f1 + f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 +
                            pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7) -
                       x * (6 * f1 + 6 * f2 - 5 * pow(pmag, 2) * pow(ppmag, 2) * f4 + 6 * pow(pmag, 2) * f5 +
                            6 * pow(ppmag, 2) * f5 + 6 * pow(pmag, 2) * f6 + 6 * pow(ppmag, 2) * f6 -
                            12 * pow(pmag, 2) * f7 + 12 * pow(ppmag, 2) * f7))) /
               2.;
    }
    if (s2.l == 2 && s1.l == 2 && s1.s == 1 && s1.j == 3)
    {
        return -0.14285714285714285 *
               (ppi *
                ((7 - 21 * pow(x, 2)) * f1 + (7 - 21 * pow(x, 2)) * f2 - 28 * pmag * ppmag * x * f3 +
                 28 * pmag * ppmag * pow(x, 3) * f3 + 5 * pow(pmag, 2) * pow(ppmag, 2) * f4 -
                 16 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 +
                 11 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 4) * f4 + pow(pmag, 2) * f5 + pow(ppmag, 2) * f5 +
                 6 * pmag * ppmag * x * f5 - 3 * pow(pmag, 2) * pow(x, 2) * f5 - 3 * pow(ppmag, 2) * pow(x, 2) * f5 -
                 10 * pmag * ppmag * pow(x, 3) * f5 + pow(pmag, 2) * f6 + pow(ppmag, 2) * f6 -
                 6 * pmag * ppmag * x * f6 - 3 * pow(pmag, 2) * pow(x, 2) * f6 - 3 * pow(ppmag, 2) * pow(x, 2) * f6 +
                 10 * pmag * ppmag * pow(x, 3) * f6 - 2 * pow(pmag, 2) * f7 + 2 * pow(ppmag, 2) * f7 +
                 6 * pow(pmag, 2) * pow(x, 2) * f7 - 6 * pow(ppmag, 2) * pow(x, 2) * f7));
    }
    if (s2.l == 2 && s1.l == 4 && s1.s == 1 && s1.j == 3)
    {
        return -(sqrt(3) * ppi *
                 (8 * pmag * ppmag * x * (-3 + 5 * pow(x, 2)) * (f5 - f6) +
                  pow(pmag, 2) * (pow(ppmag, 2) * (1 - 6 * pow(x, 2) + 5 * pow(x, 4)) * f4 +
                                  4 * (-1 + 3 * pow(x, 2)) * (f5 + f6 - 2 * f7)) +
                  pow(ppmag, 2) * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * (f5 + f6 + 2 * f7))) /
               7.;
    }
    if (s2.l == 4 && s1.l == 2 && s1.s == 1 && s1.j == 3)
    {
        return -(sqrt(3) * ppi *
                 (8 * pmag * ppmag * x * (-3 + 5 * pow(x, 2)) * (f5 - f6) +
                  pow(pmag, 2) * (pow(ppmag, 2) * (1 - 6 * pow(x, 2) + 5 * pow(x, 4)) * f4 +
                                  (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * (f5 + f6 - 2 * f7)) +
                  4 * pow(ppmag, 2) * (-1 + 3 * pow(x, 2)) * (f5 + f6 + 2 * f7))) /
               7.;
    }
    if (s2.l == 4 && s1.l == 4 && s1.s == 1 && s1.j == 3)
    {
        return (ppi *
                (7 * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * f1 + 7 * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * f2 +
                 210 * pmag * ppmag * x * f3 - 700 * pmag * ppmag * pow(x, 3) * f3 +
                 490 * pmag * ppmag * pow(x, 5) * f3 + 27 * pow(pmag, 2) * pow(ppmag, 2) * f4 -
                 267 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 +
                 485 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 4) * f4 -
                 245 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 6) * f4 - 3 * pow(pmag, 2) * f5 - 3 * pow(ppmag, 2) * f5 +
                 24 * pmag * ppmag * x * f5 + 30 * pow(pmag, 2) * pow(x, 2) * f5 + 30 * pow(ppmag, 2) * pow(x, 2) * f5 -
                 40 * pmag * ppmag * pow(x, 3) * f5 - 35 * pow(pmag, 2) * pow(x, 4) * f5 -
                 35 * pow(ppmag, 2) * pow(x, 4) * f5 - 3 * pow(pmag, 2) * f6 - 3 * pow(ppmag, 2) * f6 -
                 24 * pmag * ppmag * x * f6 + 30 * pow(pmag, 2) * pow(x, 2) * f6 + 30 * pow(ppmag, 2) * pow(x, 2) * f6 +
                 40 * pmag * ppmag * pow(x, 3) * f6 - 35 * pow(pmag, 2) * pow(x, 4) * f6 -
                 35 * pow(ppmag, 2) * pow(x, 4) * f6 + 6 * pow(pmag, 2) * f7 - 6 * pow(ppmag, 2) * f7 -
                 60 * pow(pmag, 2) * pow(x, 2) * f7 + 60 * pow(ppmag, 2) * pow(x, 2) * f7 +
                 70 * pow(pmag, 2) * pow(x, 4) * f7 - 70 * pow(ppmag, 2) * pow(x, 4) * f7)) /
               28.;
    }
    if (s2.l == 4 && s1.l == 4 && s1.s == 0 && s1.j == 4)
    {
        return (ppi * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) *
                (f1 - 3 * f2 - pow(pmag, 2) * pow(ppmag, 2) * f4 + pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                 pow(pmag, 2) * f5 - pow(ppmag, 2) * f5 - 2 * pmag * ppmag * x * f5 - pow(pmag, 2) * f6 -
                 pow(ppmag, 2) * f6 + 2 * pmag * ppmag * x * f6 + 2 * pow(pmag, 2) * f7 - 2 * pow(ppmag, 2) * f7)) /
               4.;
    }
    if (s2.l == 4 && s1.l == 4 && s1.s == 1 && s1.j == 4)
    {
        return (ppi *
                ((3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * f1 + (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * f2 +
                 6 * pmag * ppmag * x * f3 - 20 * pmag * ppmag * pow(x, 3) * f3 + 14 * pmag * ppmag * pow(x, 5) * f3 -
                 3 * pow(pmag, 2) * pow(ppmag, 2) * f4 + 27 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 2) * f4 -
                 45 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 4) * f4 +
                 21 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 6) * f4 + 3 * pow(pmag, 2) * f5 + 3 * pow(ppmag, 2) * f5 -
                 30 * pow(pmag, 2) * pow(x, 2) * f5 - 30 * pow(ppmag, 2) * pow(x, 2) * f5 -
                 40 * pmag * ppmag * pow(x, 3) * f5 + 35 * pow(pmag, 2) * pow(x, 4) * f5 +
                 35 * pow(ppmag, 2) * pow(x, 4) * f5 + 56 * pmag * ppmag * pow(x, 5) * f5 + 3 * pow(pmag, 2) * f6 +
                 3 * pow(ppmag, 2) * f6 - 30 * pow(pmag, 2) * pow(x, 2) * f6 - 30 * pow(ppmag, 2) * pow(x, 2) * f6 +
                 40 * pmag * ppmag * pow(x, 3) * f6 + 35 * pow(pmag, 2) * pow(x, 4) * f6 +
                 35 * pow(ppmag, 2) * pow(x, 4) * f6 - 56 * pmag * ppmag * pow(x, 5) * f6 - 6 * pow(pmag, 2) * f7 +
                 6 * pow(ppmag, 2) * f7 + 60 * pow(pmag, 2) * pow(x, 2) * f7 - 60 * pow(ppmag, 2) * pow(x, 2) * f7 -
                 70 * pow(pmag, 2) * pow(x, 4) * f7 + 70 * pow(ppmag, 2) * pow(x, 4) * f7)) /
               4.;
    }
    if (s2.l == 3 && s1.l == 3 && s1.s == 1 && s1.j == 4)
    {
        return -0.05555555555555555 *
               (ppi * (55 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 5) * f4 + 3 * pmag * ppmag * (9 * f3 - f5 + f6) -
                       6 * pmag * ppmag * pow(x, 2) * (27 * f3 - 5 * f5 + 5 * f6) +
                       5 * pmag * ppmag * pow(x, 4) * (27 * f3 - 7 * f5 + 7 * f6) +
                       3 * x *
                           (18 * f1 + 18 * f2 + 13 * pow(pmag, 2) * pow(ppmag, 2) * f4 + 2 * pow(pmag, 2) * f5 +
                            2 * pow(ppmag, 2) * f5 + 2 * pow(pmag, 2) * f6 + 2 * pow(ppmag, 2) * f6 -
                            4 * pow(pmag, 2) * f7 + 4 * pow(ppmag, 2) * f7) -
                       2 * pow(x, 3) *
                           (45 * f1 + 45 * f2 + 47 * pow(pmag, 2) * pow(ppmag, 2) * f4 + 5 * pow(pmag, 2) * f5 +
                            5 * pow(ppmag, 2) * f5 + 5 * pow(pmag, 2) * f6 + 5 * pow(ppmag, 2) * f6 -
                            10 * pow(pmag, 2) * f7 + 10 * pow(ppmag, 2) * f7)));
    }
    if (s2.l == 3 && s1.l == 5 && s1.s == 1 && s1.j == 4)
    {
        return -(sqrt(5) * ppi *
                 (2 * pmag * ppmag * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * (f5 - f6) +
                  pow(pmag, 2) * x *
                      (pow(ppmag, 2) * (3 - 10 * pow(x, 2) + 7 * pow(x, 4)) * f4 +
                       4 * (-3 + 5 * pow(x, 2)) * (f5 + f6 - 2 * f7)) +
                  pow(ppmag, 2) * x * (15 - 70 * pow(x, 2) + 63 * pow(x, 4)) * (f5 + f6 + 2 * f7))) /
               9.;
    }
    if (s2.l == 5 && s1.l == 3 && s1.s == 1 && s1.j == 4)
    {
        return -(sqrt(5) * ppi *
                 (2 * pmag * ppmag * (3 - 30 * pow(x, 2) + 35 * pow(x, 4)) * (f5 - f6) +
                  pow(pmag, 2) * x *
                      (pow(ppmag, 2) * (3 - 10 * pow(x, 2) + 7 * pow(x, 4)) * f4 +
                       (15 - 70 * pow(x, 2) + 63 * pow(x, 4)) * (f5 + f6 - 2 * f7)) +
                  4 * pow(ppmag, 2) * x * (-3 + 5 * pow(x, 2)) * (f5 + f6 + 2 * f7))) /
               9.;
    }
    if (s2.l == 5 && s1.l == 5 && s1.s == 1 && s1.j == 4)
    {
        return -0.027777777777777776 *
               (ppi * (-1134 * pmag * ppmag * pow(x, 6) * f3 + 567 * pow(pmag, 2) * pow(ppmag, 2) * pow(x, 7) * f4 -
                       30 * pmag * ppmag * pow(x, 2) * (27 * f3 + 2 * f5 - 2 * f6) +
                       6 * pmag * ppmag * (9 * f3 + f5 - f6) + 70 * pmag * ppmag * pow(x, 4) * (27 * f3 + f5 - f6) +
                       5 * pow(x, 3) *
                           (126 * f1 + 126 * f2 + 169 * pow(pmag, 2) * pow(ppmag, 2) * f4 - 14 * pow(pmag, 2) * f5 -
                            14 * pow(ppmag, 2) * f5 - 14 * pow(pmag, 2) * f6 - 14 * pow(ppmag, 2) * f6 +
                            28 * pow(pmag, 2) * f7 - 28 * pow(ppmag, 2) * f7) -
                       7 * pow(x, 5) *
                           (81 * f1 + 81 * f2 + 179 * pow(pmag, 2) * pow(ppmag, 2) * f4 - 9 * pow(pmag, 2) * f5 -
                            9 * pow(ppmag, 2) * f5 - 9 * pow(pmag, 2) * f6 - 9 * pow(ppmag, 2) * f6 +
                            18 * pow(pmag, 2) * f7 - 18 * pow(ppmag, 2) * f7) -
                       3 * x *
                           (45 * f1 + 45 * f2 + 53 * pow(pmag, 2) * pow(ppmag, 2) * f4 - 5 * pow(pmag, 2) * f5 -
                            5 * pow(ppmag, 2) * f5 - 5 * pow(pmag, 2) * f6 - 5 * pow(ppmag, 2) * f6 +
                            10 * pow(pmag, 2) * f7 - 10 * pow(ppmag, 2) * f7)));
    }
    else
    {
        error_message_print_abort("J out of aPWD method!");
    }
    return 0;
}

} // end namespace interaction_aPWD

#endif // INTERACTION_aPWD_HPP
