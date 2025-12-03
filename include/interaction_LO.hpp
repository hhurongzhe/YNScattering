#pragma once
#ifndef LO_INTERACTION_HPP
#define LO_INTERACTION_HPP

#include "interaction_coulomb.hpp"
#include "interaction_isospin.hpp"
#include "lib_define.hpp"
#include "state.hpp"

namespace interaction_LO
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;

// regulator function used to cut-off high momentum part in the L-S equation.
double regulator_function_order(double p1, double p2, int n, YN::YN_configs configs)
{
    double temp;
    temp = exp(-(pow(p1, n) + pow(p2, n)) / pow(configs.Lambda, n));
    return temp;
}

//------------------------------------------------------------------------------------------------------------------
// Leading order interaction in the Q = +2 particle channel:
// Sigma+ p -> Sigma+ p ,
// which is an uncoupled particle channel.
//------------------------------------------------------------------------------------------------------------------
double V_LO(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);

    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + interaction_isospin::V_isospin(4, s1, s2, configs);
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + interaction_coulomb::V_Coulomb(1, s1, s2, configs);
    }
    return regulator4 * interaction_all / twopicubic;
}
//------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------
// interaction in Q=0 coupled particle channel.
// with subscripts meaning:
// 1 : Lambda n
// 2 : Sigma0 n
// 3 : Sigma- p
//------------------------------------------------------------------------------------------------------------------

double V_LO_11(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + interaction_isospin::V_isospin(1, s1, s2, configs);
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_12(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (sqrt(1.0 / 3.0) * interaction_isospin::V_isospin(2, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_13(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (-sqrt(2.0 / 3.0) * interaction_isospin::V_isospin(2, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_21(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + sqrt(1.0 / 3.0) * interaction_isospin::V_isospin(2, s1, s2, configs);
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_22(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (1.0 / 3.0) * (interaction_isospin::V_isospin(3, s1, s2, configs) +
                                                           2.0 * interaction_isospin::V_isospin(4, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_23(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (sqrt(2.0) / 3.0) * (-interaction_isospin::V_isospin(3, s1, s2, configs) +
                                                                 interaction_isospin::V_isospin(4, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_31(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (-sqrt(2.0 / 3.0) * interaction_isospin::V_isospin(2, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_32(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (sqrt(2.0) / 3.0) * (-interaction_isospin::V_isospin(3, s1, s2, configs) +
                                                                 interaction_isospin::V_isospin(4, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + 0;
    }

    return regulator4 * interaction_all / twopicubic;
}

double V_LO_33(states::LSJ_State s1, states::LSJ_State s2, const YN::YN_configs &configs)
{
    double p1 = s1.momentum;
    double p2 = s2.momentum;
    double regulator4 = regulator_function_order(p1, p2, 4, configs);
    double interaction_all = 0.0;
    if (configs.add_strong)
    {
        interaction_all = interaction_all + (1.0 / 3.0) * (2.0 * interaction_isospin::V_isospin(3, s1, s2, configs) +
                                                           interaction_isospin::V_isospin(4, s1, s2, configs));
    }
    if (configs.add_colomb)
    {
        interaction_all = interaction_all + interaction_coulomb::V_Coulomb(-1, s1, s2, configs);
    }

    return regulator4 * interaction_all / twopicubic;
}

//------------------------------------------------------------------------------------------------------------------

} // namespace interaction_LO

#endif // LO_INTERACTION_HPP
