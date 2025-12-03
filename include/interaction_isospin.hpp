#pragma once
#ifndef INTERACTION_ISOSPIN_HPP
#define INTERACTION_ISOSPIN_HPP

#include "interaction_aPWD.hpp"
#include "lib_define.hpp"
#include "mymath.hpp"
#include "state.hpp"

namespace interaction_isospin
{
using std::string;

// Error message written and all processes stop
void error_message_print_abort(const string &error_message) { std::cout << error_message << std::endl; }

constexpr double ppi = 3.14159265358979323846;

// ome interaction in isospin space.
// isospin_channel can be 1, 2, 3, 4.
// 1: Lambda N -> Lambda N  (I=1/2)
// 2: Lambda N -> Sigma  N  (I=1/2)
// 3: Sigma  N -> Sigma  N  (I=1/2)
// 4: Sigma  N -> Sigma  N  (I=3/2)
double V_isospin_ome(int isospin_channel, states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double ff_para = configs.mb_coupling * configs.mb_coupling;
    double alpha = configs.alpha;

    // Constants for long-range meson exchange coefficients
    double long_range_pion_exchange_coeff_1 = 0.0;
    double long_range_eta_exchange_coeff_1 = (2.0 / 3.0) * (1.0 - alpha) * (4.0 * alpha - 1.0) * ff_para;
    double long_range_kaon_exchange_coeff_1 = -(1.0 / 3.0) * (1.0 + 2.0 * alpha) * (1.0 + 2.0 * alpha) * ff_para;

    double long_range_pion_exchange_coeff_2 = 2.0 * (1.0 - alpha) * ff_para;
    double long_range_eta_exchange_coeff_2 = 0.0;
    double long_range_kaon_exchange_coeff_2 = -(1.0 + 2.0 * alpha) * (1.0 - 2.0 * alpha) * ff_para;

    double long_range_pion_exchange_coeff_3 = 4.0 * alpha * ff_para;
    double long_range_eta_exchange_coeff_3 = -(2.0 / 3.0) * (1.0 - alpha) * (4.0 * alpha - 1.0) * ff_para;
    double long_range_kaon_exchange_coeff_3 = (1.0 - 2.0 * alpha) * (1.0 - 2.0 * alpha) * ff_para;

    double long_range_pion_exchange_coeff_4 = -2.0 * alpha * ff_para;
    double long_range_eta_exchange_coeff_4 = -(2.0 / 3.0) * (1.0 - alpha) * (4.0 * alpha - 1.0) * ff_para;
    double long_range_kaon_exchange_coeff_4 = -2.0 * (1.0 - 2.0 * alpha) * (1.0 - 2.0 * alpha) * ff_para;

    double pion_exchange_coeff = 0.0;
    double kaon_exchange_coeff = 0.0;
    double eta_exchange_coeff = 0.0;
    if (isospin_channel == 1)
    {
        pion_exchange_coeff = long_range_pion_exchange_coeff_1;
        kaon_exchange_coeff = long_range_kaon_exchange_coeff_1;
        eta_exchange_coeff = long_range_eta_exchange_coeff_1;
    }
    else if (isospin_channel == 2)
    {
        pion_exchange_coeff = long_range_pion_exchange_coeff_2;
        kaon_exchange_coeff = long_range_kaon_exchange_coeff_2;
        eta_exchange_coeff = long_range_eta_exchange_coeff_2;
    }
    else if (isospin_channel == 3)
    {
        pion_exchange_coeff = long_range_pion_exchange_coeff_3;
        kaon_exchange_coeff = long_range_kaon_exchange_coeff_3;
        eta_exchange_coeff = long_range_eta_exchange_coeff_3;
    }
    else if (isospin_channel == 4)
    {
        pion_exchange_coeff = long_range_pion_exchange_coeff_4;
        kaon_exchange_coeff = long_range_kaon_exchange_coeff_4;
        eta_exchange_coeff = long_range_eta_exchange_coeff_4;
    }
    else
    {
        error_message_print_abort("Not in the isospin channel table! Check the codes of interaction_isospin.hpp!\n");
    }

    std::vector<double> masslist = YN::mass_correction(isospin_channel, configs);
    double mass_pion = masslist[0];
    double mass_kaon = masslist[1];
    double mass_eta = masslist[2];

    double pmag = s1.momentum;
    double ppmag = s2.momentum;

    double temp = 0.0;
    for (int i = 0; i < 22; i = i + 1)
    {
        double x = gauss_legendre::zeros22[i];
        double w = gauss_legendre::weights22[i];
        double f1 = 0.0;
        double f2 = 0.0;
        double f3 = 0.0;
        double f4 = 0.0;
        double f5 = 0.0;
        double f6 = 0.0;
        double f7 = 0.0;
        double f8 = 0.0;

        double f6p =
            pion_exchange_coeff / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x + mass_pion * mass_pion);
        double f6e = eta_exchange_coeff / (pmag * pmag + ppmag * ppmag - 2.0 * pmag * ppmag * x + mass_eta * mass_eta);

        double f1k = kaon_exchange_coeff /
                     (pmag * pmag + ppmag * ppmag + 2.0 * pmag * ppmag * x + mass_kaon * mass_kaon) * (-0.5) *
                     (pmag * pmag + ppmag * ppmag + 2.0 * pmag * ppmag * x);
        double f2k = kaon_exchange_coeff /
                     (pmag * pmag + ppmag * ppmag + 2.0 * pmag * ppmag * x + mass_kaon * mass_kaon) * (-0.5) *
                     (-pmag * pmag - ppmag * ppmag - 2.0 * pmag * ppmag * x);
        double f5k = kaon_exchange_coeff /
                     (pmag * pmag + ppmag * ppmag + 2.0 * pmag * ppmag * x + mass_kaon * mass_kaon) * (-0.5) * (2.0);

        if (configs.exchange_pion)
        {
            f6 = f6 + f6p;
        }
        if (configs.exchange_eta)
        {
            f6 = f6 + f6e;
        }
        if (configs.exchange_kaon)
        {
            f1 = f1 + f1k;
            f2 = f2 + f2k;
            f5 = f5 + f5k;
        }

        double fa = interaction_aPWD::V_kernal_SYM(s1, s2, f1, f2, f3, f4, f5, f6, f7, f8, x);
        temp = temp + fa * w;
    }
    return temp;
}

// contact interaction in isospin space.
double V_isospin_contact(int isospin_channel, states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    if (isospin_channel == 1)
    {
        if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
        {
            return configs.C_ll_1s0; // 1S0 channel
        }
        else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
        {
            return configs.C_ll_3s1; // 3S1 channel
        }
        else
        {
            return 0.0;
        }
    }
    else if (isospin_channel == 2)
    {
        if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
        {
            return configs.C_ls_1s0; // 1S0 channel
        }
        else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
        {
            return configs.C_ls_3s1; // 3S1 channel
        }
        else
        {
            return 0.0;
        }
    }
    else if (isospin_channel == 3)
    {
        if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
        {
            return configs.CC_ss_1s0; // 1S0 channel
        }
        else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
        {
            return configs.CC_ss_3s1; // 3S1 channel
        }
        else
        {
            return 0.0;
        }
    }
    else if (isospin_channel == 4)
    {
        if (s1 == s2 && s1.l == 0 && s1.s == 0 && s1.j == 0)
        {
            return configs.C_ss_1s0; // 1S0 channel
        }
        else if (s1 == s2 && s1.l == 0 && s1.s == 1 && s1.j == 1)
        {
            return configs.C_ss_3s1; // 3S1 channel
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        return 0.0;
    }
}

// LO YN interaction in isospin space.
// chanel 1 : Lambda N -> Lambda N , I = 1/2
// chanel 2 : Lambda N -> Sigma  N , I = 1/2
// chanel 3 : Sigma  N -> Sigma  N , I = 1/2
// chanel 4 : Sigma  N -> Sigma  N , I = 3/2
double V_isospin(int isospin_channel, states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    return V_isospin_contact(isospin_channel, s1, s2, configs) + V_isospin_ome(isospin_channel, s1, s2, configs);
}

} // namespace interaction_isospin

#endif // INTERACTION_ISOSPIN_HPP
