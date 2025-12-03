#pragma once
#ifndef MATCH_HPP
#define MATCH_HPP

#include "coulomb/coulomb_calculate.hpp"
#include "lib_define.hpp"

namespace match_coulomb
{
constexpr double ppi = 3.14159265358979323846;
constexpr double twopicubic = 248.0502134423985614038105;
constexpr double fine_struct = 0.00729735253;

double match_phase_uncoupled(double deltaSL, int l, const YN::YN_configs &configs, double p_lab_now)
{
    double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
    double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
    double qq = sqrt((0.25 / ss) *
                     (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                     (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));

    double qr = qq * configs.R_Coulomb;
    double e1 = sqrt(configs.mass_sigmap * configs.mass_sigmap + qq * qq);
    double e2 = sqrt(configs.mass_proton * configs.mass_proton + qq * qq);
    double relati_corrector = (e1 * e2 + qq * qq) / (mu * (e1 + e2));
    double nn = mu / qq * fine_struct;

    double FL0, FL0p, GL0, GL0p, FL, FLp, GL, GLp, AL0, deltaCL;
    std::complex<double> eta(0, 0);
    std::complex<double> etaa(nn, 0);
    std::complex<double> z(qr, 0);

    std::vector<std::complex<double>> temp1, temp2;

    std::complex<double> l_complex(l, 0);
    temp1 = coulomb_wf_cal::coulomb_wf(l_complex, eta, z);
    temp2 = coulomb_wf_cal::coulomb_wf(l_complex, etaa, z);
    FL0 = real(temp1[0]);
    FL0p = qq * real(temp1[1]);
    GL0 = real(temp1[2]);
    GL0p = qq * real(temp1[3]);
    FL = real(temp2[0]);
    FLp = qq * real(temp2[1]);
    GL = real(temp2[2]);
    GLp = qq * real(temp2[3]);
    AL0 = (FL0 + GL0 * tan(deltaSL)) / (FL0p + GL0p * tan(deltaSL));
    deltaCL = atan((AL0 * FLp - FL) / (GL - AL0 * GLp));

    return deltaCL;
}

std::vector<double> match_phase_coupled(double rpp, double rpm, double rmm, int j, const YN::YN_configs &configs,
                                        double p_lab_now)
{
    double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
    double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
    double qq = sqrt((0.25 / ss) *
                     (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                     (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));
    double nc = mu / qq * constants::fine_struct_const;
    double qr = qq * configs.R_Coulomb;
    double factor;
    factor = -ppi * qq * mu;

    std::complex<double> eta0(0, 0);
    std::complex<double> eta(nc, 0);
    std::complex<double> z(qr, 0);

    std::vector<std::complex<double>> temp1, temp2;

    std::complex<double> lm_complex(j - 1, 0);
    temp1 = coulomb_wf_cal::coulomb_wf(lm_complex, eta0, z);
    temp2 = coulomb_wf_cal::coulomb_wf(lm_complex, eta, z);
    double F0_lm = real(temp1[0]);
    double Fp0_lm = qq * real(temp1[1]);
    double G0_lm = real(temp1[2]);
    double Gp0_lm = qq * real(temp1[3]);
    double F_lm = real(temp2[0]);
    double Fp_lm = qq * real(temp2[1]);
    double G_lm = real(temp2[2]);
    double Gp_lm = qq * real(temp2[3]);

    std::complex<double> lp_complex(j + 1, 0);
    temp1 = coulomb_wf_cal::coulomb_wf(lp_complex, eta0, z);
    temp2 = coulomb_wf_cal::coulomb_wf(lp_complex, eta, z);
    double F0_lp = real(temp1[0]);
    double Fp0_lp = qq * real(temp1[1]);
    double G0_lp = real(temp1[2]);
    double Gp0_lp = qq * real(temp1[3]);
    double F_lp = real(temp2[0]);
    double Fp_lp = qq * real(temp2[1]);
    double G_lp = real(temp2[2]);
    double Gp_lp = qq * real(temp2[3]);

    Eigen::MatrixXd F0_matrix(2, 2);
    F0_matrix << F0_lm, 0, 0, F0_lp;
    Eigen::MatrixXd Fp0_matrix(2, 2);
    Fp0_matrix << Fp0_lm, 0, 0, Fp0_lp;
    Eigen::MatrixXd G0_matrix(2, 2);
    G0_matrix << G0_lm, 0, 0, G0_lp;
    Eigen::MatrixXd Gp0_matrix(2, 2);
    Gp0_matrix << Gp0_lm, 0, 0, Gp0_lp;
    Eigen::MatrixXd F_matrix(2, 2);
    F_matrix << F_lm, 0, 0, F_lp;
    Eigen::MatrixXd Fp_matrix(2, 2);
    Fp_matrix << Fp_lm, 0, 0, Fp_lp;
    Eigen::MatrixXd G_matrix(2, 2);
    G_matrix << G_lm, 0, 0, G_lp;
    Eigen::MatrixXd Gp_matrix(2, 2);
    Gp_matrix << Gp_lm, 0, 0, Gp_lp;

    Eigen::MatrixXd RS_dimless(2, 2);
    RS_dimless << factor * rmm, factor * rpm, factor * rpm, factor * rpp;
    Eigen::MatrixXd A_matrix(2, 2);
    A_matrix = (F0_matrix + G0_matrix * RS_dimless) * ((Fp0_matrix + Gp0_matrix * RS_dimless).inverse());
    Eigen::MatrixXd Rc(2, 2);
    Rc = ((G_matrix - A_matrix * Gp_matrix).inverse()) * (A_matrix * Fp_matrix - F_matrix);
    // be careful of the possible minus sign.
    double rcpp = Rc(1, 1);
    double rcpm = Rc(0, 1);
    double rcmm = Rc(0, 0);
    double epsilon = 0.5 * atan(2.0 * rcpm / (rcmm - rcpp));
    double deltap = atan(0.5 * (rcmm + rcpp - (rcmm - rcpp) / (cos(2.0 * epsilon))));
    double deltam = atan(0.5 * (rcmm + rcpp + (rcmm - rcpp) / (cos(2.0 * epsilon))));
    double epsb = 0.5 * asin(sin(2.0 * epsilon) * sin(deltam - deltap));
    double te1 = deltap + deltam;
    double te2 = asin(tan(2.0 * epsb) / tan(2.0 * epsilon));
    double delm = 0.5 * (te1 + te2);
    double delp = 0.5 * (te1 - te2);
    std::vector<double> phase_bar_matched;
    phase_bar_matched.push_back(delm);
    phase_bar_matched.push_back(delp);
    phase_bar_matched.push_back(epsb);

    return phase_bar_matched;
}

} // end namespace match_coulomb

#endif // MATCH_HPP