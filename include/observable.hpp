#pragma once
#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include "indicators/cursor_control.hpp"
#include "indicators/progress_bar.hpp"
#include "lib_define.hpp"
#include "mmatrix.hpp"
#include "mmatrix_multichannel.hpp"
#include "phase.hpp"

namespace observable
{
using namespace indicators;

constexpr std::complex<double> im_unit{0, 1};
constexpr double ppi = 3.14159265358979323846;
constexpr double fine_struct = 0.00729735253;

void differential_cross_section(const YN::YN_configs &configs)
{
    // convert from [GeV^(-2)] to [mb], and integrating over varphi gives a 2*pi factor:
    constexpr double factor = 0.3894 * 2 * ppi;
    constexpr double angle_trans = 180.0 / ppi;
    double p_lab_now = configs.plab_dcs;
    double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
    double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
    double qq = sqrt((0.25 / ss) *
                     (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                     (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));
    double nc = mu / qq * fine_struct;

    std::vector<double> phase1, phase2;
    if (configs.match_colomb)
    {
        phase1 = phases::get_phase_uncoupled_matched(configs, p_lab_now);
        phase2 = phases::get_phase_bar_coupled_matched(configs, p_lab_now);
    }
    else
    {
        phase1 = phases::get_phase_uncoupled(configs, p_lab_now);
        phase2 = phases::get_phase_bar_coupled(configs, p_lab_now);
    }

    std::vector<double> temp;
    double costheta_min = -1;
    double costheta_max = 1;
    int interv = configs.mesh_points_integ;
    std::vector<double> costhetalist = util::make_table(costheta_min, costheta_max, interv);
    for (int i = 0; i < costhetalist.size(); i = i + 1)
    {
        double t = costhetalist[i];
        std::complex<double> M11 = mmatrix::m11(t, configs, phase1, phase2);
        std::complex<double> M10 = mmatrix::m10(t, configs, phase1, phase2);
        std::complex<double> M01 = mmatrix::m01(t, configs, phase1, phase2);
        std::complex<double> M00 = mmatrix::m00(t, configs, phase1, phase2);
        std::complex<double> Mpm = mmatrix::mpm(t, configs, phase1, phase2);
        std::complex<double> Mss = mmatrix::mss(t, configs, phase1, phase2);
        std::complex<double> Mc = -nc / (qq * (1.0 - t)) *
                                  exp(im_unit * (-nc * log(0.5 * (1.0 - t)) + 2.0 * basic_math::sigma_l_calc(0, nc)));
        double nuclear_part = 0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
        double coulomb_part = norm(Mc);
        double ci_part = (1.0 / 2.0) * real(conj(Mc) * (2.0 * M11 + M00 + Mss));
        double X = nuclear_part + coulomb_part + ci_part;
        X = X * factor;
        temp.push_back(X);
    }

    util::print_vector(costhetalist, "cos(theta)");
    util::print_vector(temp, "differential cross section (im mb)");
}

void total_cross_section(YN::YN_configs &configs)
{
    constexpr double factor = 0.3894 * 2 * ppi;
    double pstart = configs.plab_min;
    double pend = configs.plab_max;
    int interv = int((pend - pstart) / configs.plab_dis);
    double p_lab_now;

    std::vector<double> plablist = util::make_table(pstart, pend, interv);

    // Hide cursor
    show_console_cursor(false);
    // Progress Bar settings
    ProgressBar bar{option::BarWidth{50},
                    option::Start{"["},
                    option::Fill{"■"},
                    option::Lead{"■"},
                    option::Remainder{"-"},
                    option::End{" ]"},
                    option::ForegroundColor{Color::green},
                    option::ShowPercentage{true},
                    option::ShowElapsedTime{true},
                    option::ShowRemainingTime{true},
                    option::PrefixText{"Calculating total cross section "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{plablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    std::vector<double> cross_list;

    // following values are in accordence with experiment.
    double angle_max_fixed = 0.5;                // max cos(theta)
    double angle_min_fixed = -0.5;               // min cos(theta)
    int seeds_theta = configs.mesh_points_integ; // number of points in integrating over cos(theta)
    double inter = (angle_max_fixed - angle_min_fixed) / seeds_theta;

    for (int in = 0; in < plablist.size(); in = in + 1)
    {
        // bar progress update
        bar.set_progress(in);

        p_lab_now = plablist[in];
        configs.plab_dcs = p_lab_now; // update plab in configs.

        double mu = configs.mass_proton * configs.mass_sigmap / (configs.mass_proton + configs.mass_sigmap);
        double ss = configs.mass_sigmap * configs.mass_sigmap + configs.mass_proton * configs.mass_proton +
                    2.0 * configs.mass_proton * sqrt(p_lab_now * p_lab_now + configs.mass_sigmap * configs.mass_sigmap);
        double qq =
            sqrt((0.25 / ss) *
                 (ss - (configs.mass_sigmap + configs.mass_proton) * (configs.mass_sigmap + configs.mass_proton)) *
                 (ss - (configs.mass_sigmap - configs.mass_proton) * (configs.mass_sigmap - configs.mass_proton)));
        double nc = mu / qq * fine_struct;

        std::vector<double> phase1, phase2;
        if (configs.match_colomb)
        {
            phase1 = phases::get_phase_uncoupled_matched(configs, p_lab_now);
            phase2 = phases::get_phase_bar_coupled_matched(configs, p_lab_now);
        }
        else
        {
            phase1 = phases::get_phase_uncoupled(configs, p_lab_now);
            phase2 = phases::get_phase_bar_coupled(configs, p_lab_now);
        }

        double temp = 0;
        for (int ii = 0; ii < seeds_theta; ii = ii + 1)
        {
            double t = angle_min_fixed + (ii + 0.5) * inter; // cos(theta)
            std::complex<double> M11 = mmatrix::m11(t, configs, phase1, phase2);
            std::complex<double> M10 = mmatrix::m10(t, configs, phase1, phase2);
            std::complex<double> M01 = mmatrix::m01(t, configs, phase1, phase2);
            std::complex<double> M00 = mmatrix::m00(t, configs, phase1, phase2);
            std::complex<double> Mpm = mmatrix::mpm(t, configs, phase1, phase2);
            std::complex<double> Mss = mmatrix::mss(t, configs, phase1, phase2);
            std::complex<double> Mc =
                -nc / (qq * (1.0 - t)) *
                exp(im_unit * (-nc * log(0.5 * (1.0 - t)) + 2.0 * basic_math::sigma_l_calc(0, nc)));
            double nuclear_part =
                0.5 * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
            double coulomb_part = norm(Mc);
            double ci_part = (1.0 / 2.0) * real(conj(Mc) * (2.0 * M11 + M00 + Mss));
            double X = nuclear_part + coulomb_part + ci_part;
            X = X * factor;
            temp = temp + inter * X;
        }
        temp = temp * (2.0 / (angle_max_fixed - angle_min_fixed)); // normalize
        cross_list.push_back(temp);
    }
    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    std::cout.precision(6);
    util::print_vector(plablist, "p_lab (in GeV)");
    util::print_vector(cross_list, "total cross section (im mb)");
}

void differential_cross_section_multichannel(const YN::YN_configs &configs)
{
    // convert from [GeV^(-2)] to [mb], and integrating over varphi gives a 2*pi factor:
    constexpr double factor = 0.3894 * 2 * ppi;
    constexpr double angle_trans = 180.0 / ppi;
    double p_lab_now = configs.plab_dcs;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
    double qi, qf, CM_energy, ss, nc;
    switch (configs.channel_in)
    {
        case 1:
            CM_energy = sqrt(configs.mass_lambda * configs.mass_lambda + configs.mass_neutron * configs.mass_neutron +
                             2.0 * configs.mass_neutron *
                                 sqrt(configs.mass_lambda * configs.mass_lambda + p_lab_now * p_lab_now));
            ss = CM_energy * CM_energy;
            qi = sqrt(
                (0.25 / ss) *
                (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron)));
            break;
        case 2:
            CM_energy = sqrt(configs.mass_sigma0 * configs.mass_sigma0 + configs.mass_neutron * configs.mass_neutron +
                             2.0 * configs.mass_neutron *
                                 sqrt(configs.mass_sigma0 * configs.mass_sigma0 + p_lab_now * p_lab_now));
            ss = CM_energy * CM_energy;
            qi = sqrt(
                (0.25 / ss) *
                (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron)));
            break;
        case 3:
            CM_energy = sqrt(configs.mass_sigmam * configs.mass_sigmam + configs.mass_proton * configs.mass_proton +
                             2.0 * configs.mass_proton *
                                 sqrt(configs.mass_sigmam * configs.mass_sigmam + p_lab_now * p_lab_now));
            ss = CM_energy * CM_energy;
            qi = sqrt((0.25 / ss) *
                      (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                      (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton)));
            break;
        default: std::cerr << "wrong income particle channel in observable.hpp" << std::endl; std::exit(-1);
    }
    switch (configs.channel_out)
    {
        case 1:
            qf = sqrt(
                (0.25 / ss) *
                (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron)));
            break;
        case 2:
            qf = sqrt(
                (0.25 / ss) *
                (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron)));
            break;
        case 3:
            qf = sqrt((0.25 / ss) *
                      (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                      (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton)));
            break;
        default: std::cerr << "wrong outcome particle channel in observable.hpp" << std::endl; std::exit(-1);
    }
    // nc = -mu3 / qi * fine_struct * util::delta(configs.channel_in, configs.channel_out) *
    //      util::delta(configs.channel_in, 3); // coulomb effect only considered for 3 <- 3 channel.
    nc = 0;
    // double tbd = qf / qi;
    double tbd = 1;

    std::vector<std::complex<double>> Smatrix_uncoupled =
        mmatrix_multichannel::get_Smatrix_uncoupled(configs, CM_energy);
    std::vector<std::complex<double>> Smatrix_coupled = mmatrix_multichannel::get_Smatrix_coupled(configs, CM_energy);

    std::vector<double> temp;
    double costheta_min = -1;
    double costheta_max = 1;
    int interv = configs.mesh_points_integ;
    std::vector<double> costhetalist = util::make_table(costheta_min, costheta_max, interv);
    for (int i = 0; i < costhetalist.size(); i = i + 1)
    {
        double t = costhetalist[i];
        std::complex<double> M11 = mmatrix_multichannel::m11(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> M10 = mmatrix_multichannel::m10(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> M01 = mmatrix_multichannel::m01(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> M00 = mmatrix_multichannel::m00(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> Mpm = mmatrix_multichannel::mpm(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> Mss = mmatrix_multichannel::mss(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
        std::complex<double> Mc = -nc / (qi * (1.0 - t)) *
                                  exp(im_unit * (-nc * log(0.5 * (1.0 - t)) + 2.0 * basic_math::sigma_l_calc(0, nc)));
        double nuclear_part =
            0.5 * tbd * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
        double coulomb_part = norm(Mc);
        double ci_part = (1.0 / 2.0) * real(conj(Mc) * (2.0 * M11 + M00 + Mss));
        double X = nuclear_part + coulomb_part + ci_part;
        X = X * factor;
        temp.push_back(X);
    }

    std::cout.precision(3);
    std::cout << "calculating channel:           " << configs.channel_out << " <- " << configs.channel_in << "\n";
    std::cout << "laboratory momentum p_lab:     " << p_lab_now << "    GeV\n";
    std::cout << "income  c.m. momentum p_i:     " << qi << "    GeV\n";
    std::cout << "outcome c.m. momentum p_f:     " << qf << "    GeV\n\n";
    util::print_vector(costhetalist, "cos(theta)");
    std::cout.precision(6);
    util::print_vector(temp, "differential cross section (im mb)");
}

void total_cross_section_multichannel(YN::YN_configs &configs)
{
    constexpr double factor = 0.3894 * 2 * ppi;
    double mu1 = configs.mass_lambda * configs.mass_neutron / (configs.mass_lambda + configs.mass_neutron);
    double mu2 = configs.mass_sigma0 * configs.mass_neutron / (configs.mass_sigma0 + configs.mass_neutron);
    double mu3 = configs.mass_sigmam * configs.mass_proton / (configs.mass_sigmam + configs.mass_proton);
    double pstart = configs.plab_min;
    double pend = configs.plab_max;
    int interv = int((pend - pstart) / configs.plab_dis);
    double p_lab_now;

    std::vector<double> plablist = util::make_table(pstart, pend, interv);

    // Hide cursor
    show_console_cursor(false);
    // Progress Bar settings
    ProgressBar bar{option::BarWidth{50},
                    option::Start{"["},
                    option::Fill{"■"},
                    option::Lead{"■"},
                    option::Remainder{"-"},
                    option::End{" ]"},
                    option::ForegroundColor{Color::green},
                    option::ShowPercentage{true},
                    option::ShowElapsedTime{true},
                    option::ShowRemainingTime{true},
                    option::PrefixText{"Calculating total cross section "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{plablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    std::vector<double> cross_list;

    // following values are in accordence with experiment.
    double angle_max_fixed = 1.0;  // max cos(theta)
    double angle_min_fixed = -1.0; // min cos(theta)
    if (configs.channel_in == 3 && configs.channel_out == 3)
    {
        angle_max_fixed = 0.5;
        angle_min_fixed = -0.5;
    }
    int seeds_theta = configs.mesh_points_integ; // number of points in integrating over cos(theta)
    double inter = (angle_max_fixed - angle_min_fixed) / seeds_theta;

    for (int in = 0; in < plablist.size(); in = in + 1)
    {
        // bar progress update
        bar.set_progress(in);

        p_lab_now = plablist[in];
        configs.plab_dcs = p_lab_now; // update plab in configs.

        double qi, qf, CM_energy, ss, nc;
        switch (configs.channel_in)
        {
            case 1:
                CM_energy =
                    sqrt(configs.mass_lambda * configs.mass_lambda + configs.mass_neutron * configs.mass_neutron +
                         2.0 * configs.mass_neutron *
                             sqrt(configs.mass_lambda * configs.mass_lambda + p_lab_now * p_lab_now));
                ss = CM_energy * CM_energy;
                qi = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                    (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron)));
                break;
            case 2:
                CM_energy =
                    sqrt(configs.mass_sigma0 * configs.mass_sigma0 + configs.mass_neutron * configs.mass_neutron +
                         2.0 * configs.mass_neutron *
                             sqrt(configs.mass_sigma0 * configs.mass_sigma0 + p_lab_now * p_lab_now));
                ss = CM_energy * CM_energy;
                qi = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                    (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron)));
                break;
            case 3:
                CM_energy = sqrt(configs.mass_sigmam * configs.mass_sigmam + configs.mass_proton * configs.mass_proton +
                                 2.0 * configs.mass_proton *
                                     sqrt(configs.mass_sigmam * configs.mass_sigmam + p_lab_now * p_lab_now));
                ss = CM_energy * CM_energy;
                qi = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                    (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton)));
                break;
            default: std::cerr << "wrong income particle channel in observable.hpp" << std::endl; std::exit(-1);
        }
        switch (configs.channel_out)
        {
            case 1:
                qf = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_lambda + configs.mass_neutron) * (configs.mass_lambda + configs.mass_neutron)) *
                    (ss - (configs.mass_lambda - configs.mass_neutron) * (configs.mass_lambda - configs.mass_neutron)));
                break;
            case 2:
                qf = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_sigma0 + configs.mass_neutron) * (configs.mass_sigma0 + configs.mass_neutron)) *
                    (ss - (configs.mass_sigma0 - configs.mass_neutron) * (configs.mass_sigma0 - configs.mass_neutron)));
                break;
            case 3:
                qf = sqrt(
                    (0.25 / ss) *
                    (ss - (configs.mass_sigmam + configs.mass_proton) * (configs.mass_sigmam + configs.mass_proton)) *
                    (ss - (configs.mass_sigmam - configs.mass_proton) * (configs.mass_sigmam - configs.mass_proton)));
                break;
            default: std::cerr << "wrong outcome particle channel in observable.hpp" << std::endl; std::exit(-1);
        }
        // nc = -mu3 / qi * fine_struct * util::delta(configs.channel_in, configs.channel_out) *
        //      util::delta(configs.channel_in, 3); // coulomb effect only considered for 3 <- 3 channel.
        nc = 0;
        // double tbd = qf / qi;
        double tbd = 1;

        std::vector<std::complex<double>> Smatrix_uncoupled =
            mmatrix_multichannel::get_Smatrix_uncoupled(configs, CM_energy);
        std::vector<std::complex<double>> Smatrix_coupled =
            mmatrix_multichannel::get_Smatrix_coupled(configs, CM_energy);

        double temp = 0;
        for (int i = 0; i < seeds_theta; i = i + 1)
        {
            double t = angle_min_fixed + (i + 0.5) * inter; // cos(theta)
            std::complex<double> M11 =
                mmatrix_multichannel::m11(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> M10 =
                mmatrix_multichannel::m10(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> M01 =
                mmatrix_multichannel::m01(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> M00 =
                mmatrix_multichannel::m00(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> Mpm =
                mmatrix_multichannel::mpm(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> Mss =
                mmatrix_multichannel::mss(t, configs, CM_energy, Smatrix_uncoupled, Smatrix_coupled);
            std::complex<double> Mc =
                -nc / (qi * (1.0 - t)) *
                exp(im_unit * (-nc * log(0.5 * (1.0 - t)) + 2.0 * basic_math::sigma_l_calc(0, nc)));
            double nuclear_part =
                0.5 * tbd * (norm(M11) + norm(M10) + norm(Mpm) + norm(M01) + 0.5 * norm(M00) + 0.5 * norm(Mss));
            double coulomb_part = norm(Mc);
            double ci_part = (1.0 / 2.0) * real(conj(Mc) * (2.0 * M11 + M00 + Mss));
            double X = nuclear_part + coulomb_part + ci_part;
            X = X * factor;
            temp = temp + inter * X;
        }
        temp = temp * (2.0 / (angle_max_fixed - angle_min_fixed)); // normalize
        cross_list.push_back(temp);
    }
    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    std::cout << "calculating channel:           " << configs.channel_out << " <- " << configs.channel_in << "\n\n";
    std::cout.precision(6);
    util::print_vector(plablist, "p_lab (in GeV)");
    util::print_vector(cross_list, "total cross section (im mb)");
}

} // namespace observable

#endif // OBSERVABLE_HPP
