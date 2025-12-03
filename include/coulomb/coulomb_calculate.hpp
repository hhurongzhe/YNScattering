#pragma once
#ifndef COULOMB
#define COULOMB

#include "complex_functions.H"
#include "cwfcomp.cpp"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace coulomb_wf_cal
{

// regular and irregular coulomb wave function.
//[0]: F(l, eta, z)
//[1]: F'(l, eta, z)
//[2]: G(l, eta, z)
//[3]: G'(l, eta, z)
std::vector<std::complex<double>> coulomb_wf(std::complex<double> l, std::complex<double> eta, std::complex<double> z)
{
    std::vector<std::complex<double>> temp;

    bool is_renormalized = true;

    class Coulomb_wave_functions cwf(is_renormalized, l, eta);
    std::complex<double> F, dF, Hp, dHp, Hm, dHm, G, dG;
    cwf.F_dF(z, F, dF);
    cwf.G_dG(z, G, dG);

    temp.push_back(F);
    temp.push_back(dF);
    temp.push_back(G);
    temp.push_back(dG);
    return temp;
}

} // namespace coulomb_wf_cal

#endif