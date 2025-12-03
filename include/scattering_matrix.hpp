#pragma once
#ifndef SCATTERING_MATRIX_HPP
#define SCATTERING_MATRIX_HPP

#include "interaction_LO.hpp"
#include "lib_define.hpp"

namespace scattering
{
double R_uncoupled(int l, int s, int j, const YN::YN_configs &configs, double p_lab_now)
{
    double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
    double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
    double qq = sqrt((0.25 / ss) *
                     (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                     (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));
    std::vector<double> klist = configs.momentum_mesh_points;
    klist.push_back(qq);
    int N = configs.mesh_points_LS;

    // matrix equation to be solved: AX=B, (N+1)-dimension.
    Eigen::MatrixXd AA(N + 1, N + 1);
    Eigen::MatrixXd XX;
    Eigen::VectorXd BB(N + 1);

    // initialize u_j
    std::vector<double> ulist;
    double temp = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq * qq);
        double ui = 2.0 * mu * klist[i] * klist[i] * tempp;
        temp = temp + tempp;
        ulist.push_back(ui);
    }
    temp = -2.0 * mu * qq * qq * temp;
    ulist.push_back(temp);

    // initialize B vector
    states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, qq);
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[i]);
        double inter = interaction_LO::V_LO(ss1, ss2, configs);
        BB(i) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss3 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[k]);
            states::LSJ_State ss4 = states::LSJ_State::ini_LSJ_state(l, s, j, j, klist[i]);
            double inter = interaction_LO::V_LO(ss3, ss4, configs);
            double a_ik = util::delta(i, k) + ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }

    XX = AA.fullPivHouseholderQr().solve(BB);

    return XX(N, 0);
}

std::vector<double> R_coupled_triplet(int j, const YN::YN_configs &configs, double p_lab_now)
{
    double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
    double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
    double qq = sqrt((0.25 / ss) *
                     (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                     (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));
    std::vector<double> klist = configs.momentum_mesh_points;
    klist.push_back(qq);
    const int N = configs.mesh_points_LS;

    // matrix equation to be solved: AX=B, 2(N+1)-dimension.
    Eigen::MatrixXd AA(2 * (N + 1), 2 * (N + 1));
    Eigen::MatrixXd XX;
    Eigen::MatrixXd BB(2 * (N + 1), 2);

    // initialize u_j
    std::vector<double> ulist;
    double temp = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq * qq);
        double ui = 2.0 * mu * klist[i] * klist[i] * tempp;
        temp = temp + tempp;
        ulist.push_back(ui);
    }
    temp = -2.0 * mu * qq * qq * temp;
    ulist.push_back(temp);

    // initialize B matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
        double inter = interaction_LO::V_LO(ss1, ss2, configs);
        BB(i, 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i]);
        double inter = interaction_LO::V_LO(ss1, ss2, configs);
        BB(i + N + 1, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
        double inter = interaction_LO::V_LO(ss1, ss2, configs);
        BB(i, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, qq);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i]);
        ss2.momentum = klist[i];
        double inter = interaction_LO::V_LO(ss1, ss2, configs);
        BB(i + N + 1, 0) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
            double inter;
            inter = interaction_LO::V_LO(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = N + 1; i < 2 * (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i - (N + 1)]);
            double inter;
            inter = interaction_LO::V_LO(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[i]);
            double inter;
            inter = interaction_LO::V_LO(ss1, ss2, configs);
            double a_ik = ulist[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = (N + 1); i < 2 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, klist[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, klist[i - (N + 1)]);
            double inter;
            inter = interaction_LO::V_LO(ss1, ss2, configs);
            double a_ik = ulist[k] * inter;
            AA(i, k) = a_ik;
        }
    }

    XX = AA.fullPivHouseholderQr().solve(BB);

    double tpp, tpm, tmp, tmm;
    tpp = XX(N, 0);
    tpm = XX(N, 1);
    tmp = XX(2 * N + 1, 0);
    tmm = XX(2 * N + 1, 1);
    std::vector<double> tem;
    tem.push_back(tpp);
    tem.push_back(tpm);
    tem.push_back(tmp);
    tem.push_back(tmm);
    return tem;
}

} // end namespace scattering

#endif // SCATTERING_MATRIX_HPP
