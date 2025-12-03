#pragma once
#ifndef MMATRIX_MULTICHANNEL_HPP
#define MMATRIX_MULTICHANNEL_HPP

#include "S_matrix.hpp"
#include "lib_define.hpp"
#include "mmatrix.hpp"
#include "scattering_matrix_multichannel.hpp"

namespace mmatrix_multichannel
{
constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;

std::vector<std::complex<double>> get_Smatrix_uncoupled(const YN::YN_configs &configs, double CM_energy)
{
    std::vector<std::complex<double>> temp;

    int m = configs.channel_out;
    int n = configs.channel_in;
    Eigen::MatrixXd R1s0 = scattering_multichannel::R_coupledparticle_uncoupledPW(0, 0, 0, configs, CM_energy);
    Eigen::MatrixXcd S1s0 = S_matrix::get_S_matrix_uncoupledPW(R1s0, configs, CM_energy);
    std::complex<double> temp1s0 = S1s0(m - 1, n - 1);
    Eigen::MatrixXd R3p0 = scattering_multichannel::R_coupledparticle_uncoupledPW(1, 1, 0, configs, CM_energy);
    Eigen::MatrixXcd S3p0 = S_matrix::get_S_matrix_uncoupledPW(R3p0, configs, CM_energy);
    std::complex<double> temp3p0 = S3p0(m - 1, n - 1);
    temp.push_back(temp1s0);
    temp.push_back(temp3p0);

    std::complex<double> tempp0, tempp1;
    const int jmax = configs.J_uncoupled_max;
    if (jmax < 0)
    {
        std::cerr << "J_max should only be non-negative integer!" << std::endl;
        std::exit(-1);
    }
    if (jmax == 0)
    {
        // only two channels for jmax=0: 1S0, 3P0.
        return temp;
    }
    else
    {
        for (int j_temp = 1; j_temp < jmax + 1; j_temp = j_temp + 1)
        {
            Eigen::MatrixXd Rtemp0 =
                scattering_multichannel::R_coupledparticle_uncoupledPW(j_temp, 0, j_temp, configs, CM_energy);
            Eigen::MatrixXcd Stemp0 = S_matrix::get_S_matrix_uncoupledPW(Rtemp0, configs, CM_energy);
            std::complex<double> tempp0 = Stemp0(m - 1, n - 1);
            Eigen::MatrixXd Rtemp1 =
                scattering_multichannel::R_coupledparticle_uncoupledPW(j_temp, 1, j_temp, configs, CM_energy);
            Eigen::MatrixXcd Stemp1 = S_matrix::get_S_matrix_uncoupledPW(Rtemp1, configs, CM_energy);
            std::complex<double> tempp1 = Stemp1(m - 1, n - 1);
            temp.push_back(tempp0);
            temp.push_back(tempp1);
        }
    }

    return temp;
}

std::vector<std::complex<double>> get_Smatrix_coupled(const YN::YN_configs &configs, double CM_energy)
{
    std::vector<std::complex<double>> temp;

    int m = configs.channel_out;
    int n = configs.channel_in;
    const int jmax = configs.J_coupled_max;
    if (jmax < 1)
    {
        std::cerr << "error in calculation for coupled channnels!" << std::endl;
        std::exit(-1);
    }
    else
    {
        for (int j_temp = 1; j_temp < jmax + 1; j_temp = j_temp + 1)
        {
            Eigen::MatrixXd R_coupledPW =
                scattering_multichannel::R_coupledparticle_coupledPW(j_temp, configs, CM_energy);
            Eigen::MatrixXcd S_coupledPW = S_matrix::get_S_matrix_coupledPW(R_coupledPW, configs, CM_energy);
            std::complex<double> Spp = S_coupledPW(2 * (m - 1), 2 * (n - 1));
            std::complex<double> Smm = S_coupledPW(2 * (m - 1) + 1, 2 * (n - 1) + 1);
            std::complex<double> Spm = S_coupledPW(2 * (m - 1), 2 * (n - 1) + 1);
            std::complex<double> Smp = S_coupledPW(2 * (m - 1) + 1, 2 * (n - 1));
            temp.push_back(Spp);
            temp.push_back(Smm);
            temp.push_back(Spm);
            temp.push_back(Smp);
        }
    }

    return temp;
}

// scattering amplitude of "m <- n" channel.
std::complex<double> m11(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // uncoupled channels
    std::complex<double> S = Suncoupled[1];
    temp = temp + mmatrix::coeff_m11(t, 1, 1, 1, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp + 1];
        temp = temp + mmatrix::coeff_m11(t, jtemp, jtemp, 1, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    // coupled channels
    std::complex<double> Spp, Smm, Spm, Smp;
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        Spp = Scoupled[4 * jtemp - 4];
        Smm = Scoupled[4 * jtemp - 3];
        Spm = Scoupled[4 * jtemp - 2];
        Smp = Scoupled[4 * jtemp - 1];
        temp = temp +
               mmatrix::coeff_m11(t, jtemp - 1, jtemp - 1, 1, jtemp) * (Smm - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp +
               mmatrix::coeff_m11(t, jtemp + 1, jtemp + 1, 1, jtemp) * (Spp - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m11(t, jtemp - 1, jtemp + 1, 1, jtemp) * Smp / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m11(t, jtemp + 1, jtemp - 1, 1, jtemp) * Spm / (2.0 * im_unit * qq);
    }
    return temp;
}

std::complex<double> m10(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // uncoupled channels
    std::complex<double> S = Suncoupled[1];
    temp = temp + mmatrix::coeff_m10(t, 1, 1, 1, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp + 1];
        temp = temp + mmatrix::coeff_m10(t, jtemp, jtemp, 1, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    // coupled channels
    std::complex<double> Spp, Smm, Spm, Smp;
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        Spp = Scoupled[4 * jtemp - 4];
        Smm = Scoupled[4 * jtemp - 3];
        Spm = Scoupled[4 * jtemp - 2];
        Smp = Scoupled[4 * jtemp - 1];
        temp = temp +
               mmatrix::coeff_m10(t, jtemp - 1, jtemp - 1, 1, jtemp) * (Smm - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp +
               mmatrix::coeff_m10(t, jtemp + 1, jtemp + 1, 1, jtemp) * (Spp - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m10(t, jtemp - 1, jtemp + 1, 1, jtemp) * Smp / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m10(t, jtemp + 1, jtemp - 1, 1, jtemp) * Spm / (2.0 * im_unit * qq);
    }
    return temp;
}

std::complex<double> m01(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // uncoupled channels
    std::complex<double> S = Suncoupled[1];
    temp = temp + mmatrix::coeff_m01(t, 1, 1, 1, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp + 1];
        temp = temp + mmatrix::coeff_m01(t, jtemp, jtemp, 1, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    // coupled channels
    std::complex<double> Spp, Smm, Spm, Smp;
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        Spp = Scoupled[4 * jtemp - 4];
        Smm = Scoupled[4 * jtemp - 3];
        Spm = Scoupled[4 * jtemp - 2];
        Smp = Scoupled[4 * jtemp - 1];
        temp = temp +
               mmatrix::coeff_m01(t, jtemp - 1, jtemp - 1, 1, jtemp) * (Smm - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp +
               mmatrix::coeff_m01(t, jtemp + 1, jtemp + 1, 1, jtemp) * (Spp - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m01(t, jtemp - 1, jtemp + 1, 1, jtemp) * Smp / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m01(t, jtemp + 1, jtemp - 1, 1, jtemp) * Spm / (2.0 * im_unit * qq);
    }
    return temp;
}

std::complex<double> m00(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // uncoupled channels
    std::complex<double> S = Suncoupled[1];
    temp = temp + mmatrix::coeff_m00(t, 1, 1, 1, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp + 1];
        temp = temp + mmatrix::coeff_m00(t, jtemp, jtemp, 1, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    // coupled channels
    std::complex<double> Spp, Smm, Spm, Smp;
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        Spp = Scoupled[4 * jtemp - 4];
        Smm = Scoupled[4 * jtemp - 3];
        Spm = Scoupled[4 * jtemp - 2];
        Smp = Scoupled[4 * jtemp - 1];
        temp = temp +
               mmatrix::coeff_m00(t, jtemp - 1, jtemp - 1, 1, jtemp) * (Smm - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp +
               mmatrix::coeff_m00(t, jtemp + 1, jtemp + 1, 1, jtemp) * (Spp - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m00(t, jtemp - 1, jtemp + 1, 1, jtemp) * Smp / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_m00(t, jtemp + 1, jtemp - 1, 1, jtemp) * Spm / (2.0 * im_unit * qq);
    }
    return temp;
}

std::complex<double> mpm(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // uncoupled channels
    std::complex<double> S = Suncoupled[1];
    temp = temp + mmatrix::coeff_mpm(t, 1, 1, 1, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 3P0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp + 1];
        temp = temp + mmatrix::coeff_mpm(t, jtemp, jtemp, 1, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    // coupled channels
    std::complex<double> Spp, Smm, Spm, Smp;
    for (int jtemp = 1; jtemp < configs.J_coupled_max + 1; jtemp = jtemp + 1)
    {
        Spp = Scoupled[4 * jtemp - 4];
        Smm = Scoupled[4 * jtemp - 3];
        Spm = Scoupled[4 * jtemp - 2];
        Smp = Scoupled[4 * jtemp - 1];
        temp = temp +
               mmatrix::coeff_mpm(t, jtemp - 1, jtemp - 1, 1, jtemp) * (Smm - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp +
               mmatrix::coeff_mpm(t, jtemp + 1, jtemp + 1, 1, jtemp) * (Spp - util::delta(m, n)) / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_mpm(t, jtemp - 1, jtemp + 1, 1, jtemp) * Smp / (2.0 * im_unit * qq);
        temp = temp + mmatrix::coeff_mpm(t, jtemp + 1, jtemp - 1, 1, jtemp) * Spm / (2.0 * im_unit * qq);
    }
    return temp;
}

std::complex<double> mss(double t, const YN::YN_configs &configs, double CM_energy,
                         std::vector<std::complex<double>> Suncoupled, std::vector<std::complex<double>> Scoupled)
{
    int m = configs.channel_out;
    int n = configs.channel_in;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
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

    double qq;
    // switch (m * n)
    // {
    //     case 1: qq = q1; break;
    //     case 4: qq = q2; break;
    //     case 9: qq = q3; break;
    //     case 2: qq = sqrt(q1 * q2); break;
    //     case 3: qq = sqrt(q1 * q3); break;
    //     case 6: qq = sqrt(q2 * q3); break;
    //     default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    // }
    switch (n)
    {
        case 1: qq = q1; break;
        case 2: qq = q2; break;
        case 3: qq = q3; break;
        default: std::cerr << "wrong particle channel in mmatrix_multichannel.hpp" << std::endl; std::exit(-1);
    }

    std::complex<double> temp{0, 0};
    // only uncoupled channels
    std::complex<double> S = Suncoupled[0];
    temp = temp + mmatrix::coeff_mss(t, 0, 0, 0, 0) * (S - util::delta(m, n)) / (2.0 * im_unit * qq); // 1S0 channel
    for (int jtemp = 1; jtemp < configs.J_uncoupled_max + 1; jtemp = jtemp + 1)
    {
        S = Suncoupled[2 * jtemp];
        temp = temp + mmatrix::coeff_mss(t, jtemp, jtemp, 0, jtemp) * (S - util::delta(m, n)) / (2.0 * im_unit * qq);
    }
    return temp;
}

} // namespace mmatrix_multichannel

#endif // MMATRIX_MULTICHANNEL_HPP
