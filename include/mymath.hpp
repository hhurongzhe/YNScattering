#pragma once
#ifndef LO_MYMATH_HPP
#define LO_MYMATH_HPP

#include "constants.hpp"
#include "gauss_legendre/gauss_legendre.hpp"
#include <cmath>
#include <complex>
#include <vector>

namespace util
{

inline double delta(int m, int n) { return (m == n) ? 1.0 : 0.0; }

inline double sign_function(int m) { return (m % 2 == 0) ? 1.0 : -1.0; }

const std::complex<double> im_unit(0.0, 1.0); // unit imaginary number i.

std::vector<double> make_table(const double &xstart, const double &xend, const int &inter)
{
    std::vector<double> temp;
    temp.reserve(inter);
    for (int i = 0; i <= inter; i++)
    {
        double tempp = xstart + i * (xend - xstart) / inter;
        temp.push_back(tempp);
    }
    return temp;
}

void print_vector(const std::vector<double> &vec, const std::string &name)
{
    std::cout << name << ":\n";
    for (const auto &elem : vec)
    {
        std::cout << elem << ", ";
    }
    std::cout << "\n\n";
}

} // namespace util

#endif // LO_MYMATH_HPP