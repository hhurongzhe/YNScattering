#pragma once
#ifndef S_MATRIX_HPP
#define S_MATRIX_HPP

#include "lib_define.hpp"

namespace S_matrix
{
constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;

// a simple realization for kroneckerProduct of two matrix with the same dimension.
// only support complexdouble matrix.
Eigen::MatrixXcd kroneckerProduct(const Eigen::MatrixXcd &a, const Eigen::MatrixXcd &b)
{
    const int rows = a.rows();
    const int cols = a.cols();
    const int rowss = b.rows();
    const int colss = b.cols();
    if (!(rows == rowss && cols == colss))
    {
        std::cerr << "don't support kroneckerProduct of two matrix with different dimensions!" << std::endl;
        std::exit(-1);
    }

    Eigen::MatrixXcd result(rows, cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result(i, j) = a(i, j) * b(i, j);
        }
    }
    return result;
}

Eigen::MatrixXcd get_S_matrix_uncoupledPW(Eigen::MatrixXd R_uncoupled, const YN::YN_configs &configs, double CM_energy)
{
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
    double fac11 = ppi * sqrt(q1 * mu1) * sqrt(q1 * mu1);
    double fac12 = ppi * sqrt(q1 * mu1) * sqrt(q2 * mu2);
    double fac13 = ppi * sqrt(q1 * mu1) * sqrt(q3 * mu3);
    double fac21 = ppi * sqrt(q2 * mu2) * sqrt(q1 * mu1);
    double fac22 = ppi * sqrt(q2 * mu2) * sqrt(q2 * mu2);
    double fac23 = ppi * sqrt(q2 * mu2) * sqrt(q3 * mu3);
    double fac31 = ppi * sqrt(q3 * mu3) * sqrt(q1 * mu1);
    double fac32 = ppi * sqrt(q3 * mu3) * sqrt(q2 * mu2);
    double fac33 = ppi * sqrt(q3 * mu3) * sqrt(q3 * mu3);

    Eigen::MatrixXcd AA(3, 3);
    Eigen::MatrixXcd SS(3, 3);
    Eigen::MatrixXcd R_dim(3, 3);
    Eigen::MatrixXcd R_dimless(3, 3);
    Eigen::MatrixXcd TT(3, 3);
    Eigen::MatrixXcd ff(3, 3);
    auto II = Eigen::MatrixXd::Identity(3, 3);

    ff << fac11, fac12, fac13, fac21, fac22, fac23, fac31, fac32, fac33;

    R_dim = R_uncoupled;
    R_dimless = kroneckerProduct(ff, R_dim);
    AA = II + im_unit * R_dimless;
    TT = AA.inverse() * R_dimless;
    SS = II - 2.0 * im_unit * TT;

    // or we can directly get S matrix from R matrix:
    // SS = (II + im_unit * R_dimless).inverse() * (II - im_unit * R_dimless);

    return SS;
}

Eigen::MatrixXcd get_S_matrix_coupledPW(Eigen::MatrixXd R_coupled, const YN::YN_configs &configs, double CM_energy)
{
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
    double fac11 = ppi * sqrt(q1 * mu1) * sqrt(q1 * mu1);
    double fac12 = ppi * sqrt(q1 * mu1) * sqrt(q2 * mu2);
    double fac13 = ppi * sqrt(q1 * mu1) * sqrt(q3 * mu3);
    double fac21 = ppi * sqrt(q2 * mu2) * sqrt(q1 * mu1);
    double fac22 = ppi * sqrt(q2 * mu2) * sqrt(q2 * mu2);
    double fac23 = ppi * sqrt(q2 * mu2) * sqrt(q3 * mu3);
    double fac31 = ppi * sqrt(q3 * mu3) * sqrt(q1 * mu1);
    double fac32 = ppi * sqrt(q3 * mu3) * sqrt(q2 * mu2);
    double fac33 = ppi * sqrt(q3 * mu3) * sqrt(q3 * mu3);

    Eigen::MatrixXcd AA(6, 6);
    Eigen::MatrixXcd SS(6, 6);
    Eigen::MatrixXcd BB(6, 6);
    Eigen::MatrixXcd R_dim(6, 6);
    Eigen::MatrixXcd R_dimless(6, 6);
    Eigen::MatrixXcd TT(6, 6);
    auto II = Eigen::MatrixXd::Identity(6, 6);
    Eigen::MatrixXcd ff(6, 6);

    ff << fac11, fac11, fac12, fac12, fac13, fac13, fac11, fac11, fac12, fac12, fac13, fac13, fac21, fac21, fac22,
        fac22, fac23, fac23, fac21, fac21, fac22, fac22, fac23, fac23, fac31, fac31, fac32, fac32, fac33, fac33, fac31,
        fac31, fac32, fac32, fac33, fac33;

    R_dim = R_coupled;
    R_dimless = kroneckerProduct(ff, R_dim);
    AA = II + im_unit * R_dimless;
    TT = AA.inverse() * R_dimless;
    SS = II - 2.0 * im_unit * TT;

    // or we can directly get S matrix from R matrix:
    // SS = (II + im_unit * R_dimless).inverse() * (II - im_unit * R_dimless);

    return SS;
}
} // end namespace S_matrix

#endif // S_MATRIX_HPP