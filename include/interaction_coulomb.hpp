#pragma once
#ifndef COULOMB_INTERACTION_HPP
#define COULOMB_INTERACTION_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "mymath.hpp"
#include "state.hpp"

namespace interaction_coulomb
{

// function W(t) in Coulomb interaction in momentum space.
double Coulomb_wt(double p1, double p2, double tt, double rr)
{
    double temp;
    double tem = p1 * p1 + p2 * p2 - 2.0 * p1 * p2 * tt;
    temp = (1.0 - cos(sqrt(tem) * rr)) / tem;
    return temp;
}

// cut-off colomb interaction in momentum space,
// partial-wave matrix elements in LSJ-SYM basis.
double V_Coulomb(double charge_factor, states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double temp = 0.0;

    double parameter_colomb = charge_factor * configs.C_C;

    double pmag = s1.momentum;
    double ppmag = s2.momentum;

    double temp_fac = pmag * pmag / configs.mass_proton / configs.mass_proton;
    double relati_corrector = (1.0 + 2.0 * temp_fac) / (sqrt(1.0 + temp_fac));

    for (int i = 0; i < 22; i = i + 1)
    {
        double x = gauss_legendre::zeros22[i];
        double w = gauss_legendre::weights22[i];
        double f1_colomb = parameter_colomb * Coulomb_wt(pmag, ppmag, x, configs.R_Coulomb);

        double fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1_colomb, 0, 0, 0, 0, 0, 0, 0, x);
        temp = temp + fa * w;
    }
    return temp;
}

} // namespace interaction_coulomb

#endif // COULOMB_INTERACTION_HPP
