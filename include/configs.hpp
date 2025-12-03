#pragma once
#ifndef YN_CONFIGS_HPP
#define YN_CONFIGS_HPP

#include "mymath.hpp"
#include "util/inifile.hpp"

namespace YN
{

struct YN_configs
{

    // ***** interaction section *****
    double axial_current_coupling_constant;
    double pion_decay_constant;
    double alpha;
    double mb_coupling;
    double C_ll_1s0, C_ll_3s1, C_ss_1s0, C_ss_3s1, C_ls_1s0, C_ls_3s1, CC_ss_1s0, CC_ss_3s1;
    double Lambda;
    bool exchange_pion;
    bool exchange_kaon;
    bool exchange_eta;
    bool loop_corr;
    bool add_strong;
    bool add_colomb;
    bool match_colomb;

    // ***** kinetic parameters section ****
    double plab_min;
    double plab_max;
    double plab_dis;
    double plab_dcs;
    int channel_in;
    int channel_out;

    // ***** meson masses section ****
    double mass_pion_charged;
    double mass_pion_neutral;
    double mass_pion_averaged;
    double mass_kaon;
    double mass_eta;

    // ***** baryon masses section ****
    double mass_proton;
    double mass_neutron;
    double mass_nucleon;
    double mass_sigmap;
    double mass_sigmam;
    double mass_sigma0;
    double mass_sigma_averaged;
    double mass_lambda;

    // ***** numerical parameters section ****
    int mesh_points_LS;
    int mesh_points_integ;
    int mesh_points_interaction;
    std::vector<double> momentum_mesh_points;
    std::vector<double> momentum_mesh_weights;

    double C_Gaussian;
    double R_Coulomb;
    double C_C; // constant C_C of coulomb force in momentum space.

    int J_uncoupled_max;
    int J_coupled_max;

    // ***** output section *****
    std::string result_dir;
    std::string result_name;

    // constructor
    YN_configs(const inifile_system::inifile &ini);

    // generate a result file name
    std::string result_file() const;
};

YN_configs::YN_configs(const inifile_system::inifile &ini)
{

    // ***** interaction section *****
    auto sec = ini.section("interaction");
    axial_current_coupling_constant = sec.get_double("axial_current_coupling_constant");
    pion_decay_constant = 0.001 * sec.get_double("pion_decay_constant");
    alpha = sec.get_double("alpha");
    mb_coupling = axial_current_coupling_constant / (2.0 * pion_decay_constant);

    C_ll_1s0 = 10000.0 * sec.get_double("C_ll_1s0");
    C_ll_3s1 = 10000.0 * sec.get_double("C_ll_3s1");
    C_ss_1s0 = 10000.0 * sec.get_double("C_ss_1s0");
    C_ss_3s1 = 10000.0 * sec.get_double("C_ss_3s1");
    C_ls_3s1 = 10000.0 * sec.get_double("C_ls_3s1");
    C_ls_1s0 = 3.0 * (C_ll_1s0 - C_ss_1s0);
    CC_ss_1s0 = 9.0 * C_ll_1s0 - 8.0 * C_ss_1s0;
    CC_ss_3s1 = C_ll_3s1;

    Lambda = 0.001 * sec.get_double("Lambda");

    exchange_pion = sec.get_bool("exchange_pion");
    exchange_kaon = sec.get_bool("exchange_kaon");
    exchange_eta = sec.get_bool("exchange_eta");
    loop_corr = sec.get_bool("loop_corr");
    add_strong = sec.get_bool("add_strong");
    add_colomb = sec.get_bool("add_colomb");
    match_colomb = sec.get_bool("match_colomb");

    // ***** kinetic parameters section ****
    sec = ini.section("kinetics");
    plab_min = 0.001 * sec.get_double("plab_min");
    plab_max = 0.001 * sec.get_double("plab_max");
    plab_dis = 0.001 * sec.get_double("plab_dis");
    plab_dcs = 0.001 * sec.get_double("plab_dcs");
    channel_in = sec.get_int("channel_in");
    channel_out = sec.get_int("channel_out");

    // ***** meson masses section ****
    sec = ini.section("meson-masses");
    mass_pion_charged = 0.001 * sec.get_double("mass_pion_charged");
    mass_pion_neutral = 0.001 * sec.get_double("mass_pion_neutral");
    mass_pion_averaged = 0.001 * sec.get_double("mass_pion_averaged");
    mass_kaon = 0.001 * sec.get_double("mass_kaon");
    mass_eta = 0.001 * sec.get_double("mass_eta");

    // ***** baryon masses section ****
    sec = ini.section("baryon-masses");
    mass_proton = 0.001 * sec.get_double("mass_proton");
    mass_neutron = 0.001 * sec.get_double("mass_neutron");
    mass_nucleon = 0.001 * sec.get_double("mass_nucleon");
    mass_sigmap = 0.001 * sec.get_double("mass_sigmap");
    mass_sigmam = 0.001 * sec.get_double("mass_sigmam");
    mass_sigma0 = 0.001 * sec.get_double("mass_sigma0");
    mass_sigma_averaged = 0.001 * sec.get_double("mass_sigma_averaged");
    mass_lambda = 0.001 * sec.get_double("mass_lambda");

    // ***** numerical parameters section ****
    sec = ini.section("numerical-parameters");
    mesh_points_LS = sec.get_int("mesh_points_LS");
    mesh_points_integ = sec.get_int("mesh_points_integ");
    mesh_points_interaction = sec.get_int("mesh_points_interaction");
    C_Gaussian = 0.001 * sec.get_double("C_Gaussian");
    R_Coulomb = 5.06773125 * sec.get_double("R_Coulomb");
    C_C = 0.0917012;
    J_uncoupled_max = sec.get_int("J_uncoupled_max");
    J_coupled_max = sec.get_int("J_coupled_max");
    momentum_mesh_points = gauss_legendre::make_klist(mesh_points_LS, C_Gaussian);
    momentum_mesh_weights = gauss_legendre::make_slist(mesh_points_LS, C_Gaussian);

    // ***** output section *****
    sec = ini.section("output");
    result_dir = sec.get_string("result_dir");
    result_name = sec.get_string("result_name");
    if (result_dir.back() != '/')
    {
        result_dir += "/";
    }
};

std::string file_stem(const std::string &file)
{
    auto p1 = file.find_last_of('/') + 1;
    auto p2 = file.find_last_of('.');
    return file.substr(p1, p2 - p1);
}

std::string YN_configs::result_file() const
{
    std::ostringstream oss;
    oss << result_dir << result_name << ".trace";
    return oss.str();
}

std::vector<double> mass_correction(int isospin_channel, const YN::YN_configs &configs)
{
    std::vector<double> temp;
    double m_pi = configs.mass_pion_averaged;
    double m_lambda = configs.mass_lambda;
    double m_nucleon = configs.mass_nucleon;
    double m_k = configs.mass_kaon;
    double m_eta = configs.mass_eta;
    double m_sigma = configs.mass_sigma_averaged;
    bool corr = configs.loop_corr;
    double mp, mk, me;
    if (corr)
    {
        if (isospin_channel == 1)
        {
            mp = m_pi;
            double dm2 = (m_lambda - m_nucleon) * (m_lambda - m_nucleon);
            mk = sqrt(m_k * m_k - dm2);
            me = m_eta;
        }
        else if (isospin_channel == 2)
        {
            double dmp2 = (m_lambda - m_sigma) * (m_lambda - m_sigma) / 4.0;
            mp = sqrt(m_pi * m_pi - dmp2);
            double dmk2 = ((m_lambda + m_sigma) / 2.0 - m_nucleon) * ((m_lambda + m_sigma) / 2.0 - m_nucleon);
            mk = sqrt(m_k * m_k - dmk2);
            me = m_eta;
        }
        else if (isospin_channel == 3)
        {
            mp = m_pi;
            double dm2 = (m_sigma - m_nucleon) * (m_sigma - m_nucleon);
            mk = sqrt(m_k * m_k - dm2);
            me = m_eta;
        }
        else if (isospin_channel == 4)
        {
            mp = m_pi;
            double dm2 = (m_sigma - m_nucleon) * (m_sigma - m_nucleon);
            mk = sqrt(m_k * m_k - dm2);
            me = m_eta;
        }
        else
        {
            std::cout << "Not in the isospin channel table!\n Check the codes!\n" << std::endl;
        }
    }
    else
    {
        mp = m_pi;
        mk = m_k;
        me = m_eta;
    }
    temp.push_back(mp);
    temp.push_back(mk);
    temp.push_back(me);
    return temp;
}

} // namespace YN

#endif // YN_CONFIGS_HPP