#pragma once
#ifndef SCATTERING_MATRIX_MULTICHANNEL_HPP
#define SCATTERING_MATRIX_MULTICHANNEL_HPP

#include "interaction_LO.hpp"
#include "lib_define.hpp"

namespace scattering_multichannel
{

// solving coupled-channel LS equation for coupled particle channel.
// particle channel index are:
// 1: Lambda-Neutron
// 2: Sigma0-Neutron
// 3: Sigmam-Protron
Eigen::MatrixXd R_coupledparticle_uncoupledPW(int l, int s, int j, const YN::YN_configs &configs, double CM_energy)
{
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
    int N = configs.mesh_points_LS;
    std::vector<double> klist = configs.momentum_mesh_points;

    double ss = CM_energy * CM_energy;
    double qq1 = (0.25 / ss) *
                 (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                 (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron));
    double qq2 = (0.25 / ss) *
                 (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                 (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron));
    double qq3 = (0.25 / ss) *
                 (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                 (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton));
    double q1, q2, q3;
    q1 = sqrt(qq1);
    q2 = sqrt(qq2);
    q3 = sqrt(qq3);

    std::vector<double> k1list, k2list, k3list;
    for (int i = 0; i < N; i = i + 1)
    {
        k1list.push_back(klist[i]);
        k2list.push_back(klist[i]);
        k3list.push_back(klist[i]);
    }
    k1list.push_back(q1);
    k2list.push_back(q2);
    k3list.push_back(q3);

    // matrix equation to be solved: AX=B, 3(N+1)-dimension.
    Eigen::MatrixXd AA(3 * (N + 1), 3 * (N + 1));
    Eigen::MatrixXd XX;
    Eigen::MatrixXd BB(3 * (N + 1), 3);

    // initialize u_j
    std::vector<double> ulist1, ulist2, ulist3;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp1 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq1);
        double tempp2 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq2);
        double tempp3 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq3);
        double ui1 = 2.0 * mu1 * klist[i] * klist[i] * tempp1;
        double ui2 = 2.0 * mu2 * klist[i] * klist[i] * tempp2;
        double ui3 = 2.0 * mu3 * klist[i] * klist[i] * tempp3;
        temp1 = temp1 + tempp1;
        temp2 = temp2 + tempp2;
        temp3 = temp3 + tempp3;
        ulist1.push_back(ui1);
        ulist2.push_back(ui2);
        ulist3.push_back(ui3);
    }
    ulist1.push_back(-2.0 * mu1 * qq1 * temp1);
    ulist2.push_back(-2.0 * mu2 * qq2 * temp2);
    ulist3.push_back(-2.0 * mu3 * qq3 * temp3);

    // initialize B matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

        double inter = interaction_LO::V_LO_11(ss1, ss2, configs);
        BB(i, 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i]);

        double inter = interaction_LO::V_LO_21(ss1, ss2, configs);
        BB(i + N + 1, 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i]);

        double inter = interaction_LO::V_LO_31(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

        double inter = interaction_LO::V_LO_12(ss1, ss2, configs);
        BB(i, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i]);

        double inter = interaction_LO::V_LO_22(ss1, ss2, configs);
        BB(i + N + 1, 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i]);

        double inter = interaction_LO::V_LO_32(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 1) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

        double inter = interaction_LO::V_LO_13(ss1, ss2, configs);
        BB(i, 2) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i]);

        double inter = interaction_LO::V_LO_23(ss1, ss2, configs);
        BB(i + N + 1, 2) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i]);

        double inter = interaction_LO::V_LO_33(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 2) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

            double inter = interaction_LO::V_LO_11(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist1[k] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = N + 1; i < 2 * (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i - (N + 1)]);

            double inter = interaction_LO::V_LO_22(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist2[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i - 2 * (N + 1)]);

            double inter = interaction_LO::V_LO_33(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist3[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = N + 1; k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

            double inter = interaction_LO::V_LO_12(ss1, ss2, configs);
            double a_ik = ulist2[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[i]);

            double inter = interaction_LO::V_LO_13(ss1, ss2, configs);
            double a_ik = ulist3[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = (N + 1); i < 2 * (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i - (N + 1)]);

            double inter = interaction_LO::V_LO_23(ss1, ss2, configs);
            double a_ik = ulist3[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = (N + 1); i < 2 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[i - (N + 1)]);

            double inter = interaction_LO::V_LO_21(ss1, ss2, configs);
            double a_ik = ulist1[k] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i - 2 * (N + 1)]);

            double inter = interaction_LO::V_LO_31(ss1, ss2, configs);
            double a_ik = ulist1[k] * inter;
            AA(i, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = (N + 1); k < 2 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k2list[k - (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(l, s, j, j, k3list[i - 2 * (N + 1)]);

            double inter = interaction_LO::V_LO_32(ss1, ss2, configs);
            double a_ik = ulist2[k - (N + 1)] * inter;
            AA(i, k) = a_ik;
        }
    }

    XX = AA.fullPivHouseholderQr().solve(BB);

    Eigen::MatrixXd tem(3, 3);
    double r11, r12, r13, r21, r22, r23, r31, r32, r33;
    r11 = XX(N, 0);
    r12 = XX(N, 1);
    r13 = XX(N, 2);
    r21 = XX(2 * N + 1, 0);
    r22 = XX(2 * N + 1, 1);
    r23 = XX(2 * N + 1, 2);
    r31 = XX(3 * N + 2, 0);
    r32 = XX(3 * N + 2, 1);
    r33 = XX(3 * N + 2, 2);
    tem << r11, r12, r13, r21, r22, r23, r31, r32, r33;
    return tem;
}

Eigen::MatrixXd R_coupledparticle_coupledPW(int j, const YN::YN_configs &configs, double CM_energy)
{
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
    int N = configs.mesh_points_LS;
    std::vector<double> klist = configs.momentum_mesh_points;

    double ss = CM_energy * CM_energy;
    double qq1 = (0.25 / ss) *
                 (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                 (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron));
    double qq2 = (0.25 / ss) *
                 (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                 (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron));
    double qq3 = (0.25 / ss) *
                 (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                 (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton));
    double q1, q2, q3;
    q1 = sqrt(qq1);
    q2 = sqrt(qq2);
    q3 = sqrt(qq3);

    std::vector<double> k1list, k2list, k3list;
    for (int i = 0; i < N; i = i + 1)
    {
        k1list.push_back(klist[i]);
        k2list.push_back(klist[i]);
        k3list.push_back(klist[i]);
    }
    k1list.push_back(q1);
    k2list.push_back(q2);
    k3list.push_back(q3);

    // matrix equation to be solved: AX=B, 6(N+1)-dimension.
    Eigen::MatrixXd AA(6 * (N + 1), 6 * (N + 1));
    Eigen::MatrixXd XX;
    Eigen::MatrixXd BB(6 * (N + 1), 6);

    // initialize u_j
    std::vector<double> ulist1, ulist2, ulist3;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    for (int i = 0; i < N; i = i + 1)
    {
        double tempp1 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq1);
        double tempp2 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq2);
        double tempp3 = configs.momentum_mesh_weights[i] / (klist[i] * klist[i] - qq3);
        double ui1 = 2.0 * mu1 * klist[i] * klist[i] * tempp1;
        double ui2 = 2.0 * mu2 * klist[i] * klist[i] * tempp2;
        double ui3 = 2.0 * mu3 * klist[i] * klist[i] * tempp3;
        temp1 = temp1 + tempp1;
        temp2 = temp2 + tempp2;
        temp3 = temp3 + tempp3;
        ulist1.push_back(ui1);
        ulist2.push_back(ui2);
        ulist3.push_back(ui3);
    }
    ulist1.push_back(-2.0 * mu1 * qq1 * temp1);
    ulist2.push_back(-2.0 * mu2 * qq2 * temp2);
    ulist3.push_back(-2.0 * mu3 * qq3 * temp3);

    // initialize B matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        double inter = interaction_LO::V_LO_11(ss1, ss2, configs);
        BB(i, 0) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_11(ss1, ss2, configs);
        BB(i + N + 1, 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_11(ss1, ss2, configs);
        BB(i, 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_11(ss1, ss2, configs);
        BB(i + N + 1, 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        double inter = interaction_LO::V_LO_21(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 0) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_21(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_21(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_21(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        double inter = interaction_LO::V_LO_31(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 0) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_31(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_31(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 1) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q1);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_31(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 0) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        double inter = interaction_LO::V_LO_12(ss1, ss2, configs);
        BB(i, 2) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_12(ss1, ss2, configs);
        BB(i + N + 1, 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_12(ss1, ss2, configs);
        BB(i, 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_12(ss1, ss2, configs);
        BB(i + N + 1, 2) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        double inter = interaction_LO::V_LO_22(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 2) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_22(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_22(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_22(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 2) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        double inter = interaction_LO::V_LO_32(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 2) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_32(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_32(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 3) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q2);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_32(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 2) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        double inter = interaction_LO::V_LO_13(ss1, ss2, configs);
        BB(i, 4) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_13(ss1, ss2, configs);
        BB(i + N + 1, 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_13(ss1, ss2, configs);
        BB(i, 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
        inter = interaction_LO::V_LO_13(ss1, ss2, configs);
        BB(i + N + 1, 4) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        double inter = interaction_LO::V_LO_23(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 4) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_23(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_23(ss1, ss2, configs);
        BB(i + 2 * (N + 1), 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i]);
        inter = interaction_LO::V_LO_23(ss1, ss2, configs);
        BB(i + 3 * (N + 1), 4) = inter;
    }
    for (int i = 0; i < N + 1; i = i + 1)
    {
        states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        double inter = interaction_LO::V_LO_33(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 4) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_33(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_33(ss1, ss2, configs);
        BB(i + 4 * (N + 1), 5) = inter;

        ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, q3);
        ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i]);
        inter = interaction_LO::V_LO_33(ss1, ss2, configs);
        BB(i + 5 * (N + 1), 4) = inter;
    }

    // initialize A matrix
    for (int i = 0; i < N + 1; i = i + 1)
    {
        for (int k = 0; k < N + 1; k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            double inter = interaction_LO::V_LO_11(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist1[k] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_11(ss1, ss2, configs);
            a_ik = util::delta(i, k) + ulist1[k] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_11(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_11(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            double inter = interaction_LO::V_LO_22(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_22(ss1, ss2, configs);
            a_ik = util::delta(i, k) + ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_22(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_22(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 4 * (N + 1); i < 5 * (N + 1); i = i + 1)
    {
        for (int k = 4 * (N + 1); k < 5 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            double inter = interaction_LO::V_LO_33(ss1, ss2, configs);
            double a_ik = util::delta(i, k) + ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_33(ss1, ss2, configs);
            a_ik = util::delta(i, k) + ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_33(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_33(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            double inter = interaction_LO::V_LO_12(ss1, ss2, configs);
            double a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_12(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_12(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_12(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 0; i < (N + 1); i = i + 1)
    {
        for (int k = 4 * (N + 1); k < 5 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            double inter = interaction_LO::V_LO_13(ss1, ss2, configs);
            double a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_13(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_13(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[i]);
            inter = interaction_LO::V_LO_13(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = 4 * (N + 1); k < 5 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            double inter = interaction_LO::V_LO_23(ss1, ss2, configs);
            double a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_23(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_23(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[k - 4 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_23(ss1, ss2, configs);
            a_ik = ulist3[k - 4 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 2 * (N + 1); i < 3 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            double inter = interaction_LO::V_LO_21(ss1, ss2, configs);
            double a_ik = ulist1[k] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_21(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_21(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[i - 2 * (N + 1)]);
            inter = interaction_LO::V_LO_21(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 4 * (N + 1); i < 5 * (N + 1); i = i + 1)
    {
        for (int k = 0; k < (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            double inter = interaction_LO::V_LO_31(ss1, ss2, configs);
            double a_ik = ulist1[k] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_31(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_31(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k1list[k]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_31(ss1, ss2, configs);
            a_ik = ulist1[k] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }
    for (int i = 4 * (N + 1); i < 5 * (N + 1); i = i + 1)
    {
        for (int k = 2 * (N + 1); k < 3 * (N + 1); k = k + 1)
        {
            states::LSJ_State ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            states::LSJ_State ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            double inter = interaction_LO::V_LO_32(ss1, ss2, configs);
            double a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_32(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_32(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i, k + N + 1) = a_ik;

            ss1 = states::LSJ_State::ini_LSJ_state(j + 1, 1, j, j, k2list[k - 2 * (N + 1)]);
            ss2 = states::LSJ_State::ini_LSJ_state(j - 1, 1, j, j, k3list[i - 4 * (N + 1)]);
            inter = interaction_LO::V_LO_32(ss1, ss2, configs);
            a_ik = ulist2[k - 2 * (N + 1)] * inter;
            AA(i + N + 1, k) = a_ik;
        }
    }

    XX = AA.fullPivHouseholderQr().solve(BB);

    double t11pp, t11pm, t11mp, t11mm, t12pp, t12pm, t12mp, t12mm, t13pp, t13pm, t13mp, t13mm;
    double t21pp, t21pm, t21mp, t21mm, t22pp, t22pm, t22mp, t22mm, t23pp, t23pm, t23mp, t23mm;
    double t31pp, t31pm, t31mp, t31mm, t32pp, t32pm, t32mp, t32mm, t33pp, t33pm, t33mp, t33mm;
    t11pp = XX(N, 0);
    t11pm = XX(N, 1);
    t11mp = XX(N + (N + 1), 0);
    t11mm = XX(N + (N + 1), 1);
    t21pp = XX(N + (N + 1), 0);
    t21pm = XX(N + (N + 1), 1);
    t21mp = XX(N + 2 * (N + 1), 0);
    t21mm = XX(N + 2 * (N + 1), 1);
    t31pp = XX(N + 3 * (N + 1), 0);
    t31pm = XX(N + 3 * (N + 1), 1);
    t31mp = XX(N + 4 * (N + 1), 0);
    t31mm = XX(N + 4 * (N + 1), 1);
    t12pp = XX(N, 2);
    t12pm = XX(N, 3);
    t12mp = XX(N + (N + 1), 2);
    t12mm = XX(N + (N + 1), 3);
    t22pp = XX(N + (N + 1), 2);
    t22pm = XX(N + (N + 1), 3);
    t22mp = XX(N + 2 * (N + 1), 2);
    t22mm = XX(N + 2 * (N + 1), 3);
    t32pp = XX(N + 3 * (N + 1), 2);
    t32pm = XX(N + 3 * (N + 1), 3);
    t32mp = XX(N + 4 * (N + 1), 2);
    t32mm = XX(N + 4 * (N + 1), 3);
    t13pp = XX(N, 4);
    t13pm = XX(N, 5);
    t13mp = XX(N + (N + 1), 4);
    t13mm = XX(N + (N + 1), 5);
    t23pp = XX(N + (N + 1), 4);
    t23pm = XX(N + (N + 1), 5);
    t23mp = XX(N + 2 * (N + 1), 4);
    t23mm = XX(N + 2 * (N + 1), 5);
    t33pp = XX(N + 3 * (N + 1), 4);
    t33pm = XX(N + 3 * (N + 1), 5);
    t33mp = XX(N + 4 * (N + 1), 4);
    t33mm = XX(N + 4 * (N + 1), 5);

    Eigen::MatrixXd tem(6, 6);
    tem << t11pp, t11pm, t12pp, t12pm, t13pp, t13pm, t11mp, t11mm, t12mp, t12mm, t13mp, t13mm, t21pp, t21pm, t22pp,
        t22pm, t23pp, t23pm, t21mp, t21mm, t22mp, t22mm, t23mp, t23mm, t31pp, t31pm, t32pp, t32pm, t33pp, t33pm, t31mp,
        t31mm, t32mp, t32mm, t33mp, t33mm;

    return tem;
}

} // namespace scattering_multichannel

#endif // SCATTERING_MATRIX_MULTICHANNEL_HPP
