#pragma once
#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include "indicators/cursor_control.hpp"
#include "indicators/progress_bar.hpp"
#include "lib_define.hpp"
#include "phase.hpp"

namespace evaluate
{
using namespace indicators;

// phase shifts evaluator for scattering:
// Sigma+ + p -> Sigma+ + p ,
// which not considering coulomb interaction.
void evaluate_phase_unmatched(const YN::YN_configs &configs)
{
    double pstart = configs.plab_min;
    double pend = configs.plab_max;
    int interv = int((pend - pstart) / configs.plab_dis);
    int de = configs.mesh_points_LS;
    double angle_trans = 180.0 / constants::pi;
    double p_lab_now;

    std::vector<double> plablist = util::make_table(pstart, pend, interv);

    //* store calculated phase shifts in matrix by Eigen:
    Eigen::MatrixXd phase_uncoupled(plablist.size(), 2 * configs.J_uncoupled_max + 2);
    Eigen::MatrixXd phase_coupled(plablist.size(), 3 * configs.J_uncoupled_max);

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
                    option::PrefixText{"Calculating phase shifts "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{plablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    // note that in, ii and  jj are dummy variables.
    for (int in = 0; in < plablist.size(); in = in + 1)
    {
        p_lab_now = plablist[in];

        // bar progress update
        bar.set_progress(in);

        std::vector<double> phase1 = phases::get_phase_uncoupled(configs, p_lab_now);
        std::vector<double> phase2 = phases::get_phase_bar_coupled(configs, p_lab_now);
        for (int ii = 0; ii < 2 * configs.J_uncoupled_max + 2; ii = ii + 1)
        {
            phase_uncoupled(in, ii) = angle_trans * phase1[ii];
        }
        for (int ii = 0; ii < 3 * configs.J_coupled_max; ii = ii + 1)
        {
            phase_coupled(in, ii) = angle_trans * phase2[ii];
        }
    }

    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    // write results in the output file.
    std::ofstream fp(configs.result_file());
    std::cout.precision(6);
    auto name_uncoupled = phases::partial_wave_name_uncoupled(configs);
    auto name_coupled = phases::partial_wave_name_coupled(configs);
    // print partial wave names on the first line in output file:
    fp << std::fixed << "p_lab";
    for (const auto &name : name_uncoupled)
    {
        fp << "\t" << std::setw(15) << name;
    }
    for (const auto &name : name_coupled)
    {
        fp << "\t" << std::setw(15) << name;
    }
    fp << "\n";
    for (int ii = 0; ii < plablist.size(); ii = ii + 1)
    {
        fp << std::fixed << 1000.0 * plablist[ii];
        for (int jj = 0; jj < 2 * configs.J_uncoupled_max + 2; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_uncoupled(ii, jj);
        }
        for (int jj = 0; jj < 3 * configs.J_coupled_max; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_coupled(ii, jj);
        }
        fp << '\n';
    }
    fp.close();
}

// phase shifts evaluator for scattering:
// Sigma+ + p -> Sigma+ + p ,
// which fully considering coulomb interaction.
void evaluate_phase_matched(const YN::YN_configs &configs)
{
    double pstart = configs.plab_min;
    double pend = configs.plab_max;
    int interv = int((pend - pstart) / configs.plab_dis);
    int de = configs.mesh_points_LS;
    double angle_trans = 180.0 / constants::pi;
    double p_lab_now;

    std::vector<double> plablist = util::make_table(pstart, pend, interv);

    //* store calculated phase shifts in matrix by Eigen:
    Eigen::MatrixXd phase_uncoupled(plablist.size(), 2 * configs.J_uncoupled_max + 2);
    Eigen::MatrixXd phase_coupled(plablist.size(), 3 * configs.J_uncoupled_max);

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
                    option::PrefixText{"Calculating phase shifts "},
                    option::ProgressType{ProgressType::incremental},
                    option::MaxProgress{plablist.size()},
                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

    // note that in, ii and  jj are dummy variables.
    for (int in = 0; in < plablist.size(); in = in + 1)
    {
        p_lab_now = plablist[in];

        // bar progress update
        bar.set_progress(in);

        std::vector<double> phase1 = phases::get_phase_uncoupled_matched(configs, p_lab_now);
        std::vector<double> phase2 = phases::get_phase_bar_coupled_matched(configs, p_lab_now);
        for (int ii = 0; ii < 2 * configs.J_uncoupled_max + 2; ii = ii + 1)
        {
            phase_uncoupled(in, ii) = angle_trans * phase1[ii];
        }
        for (int ii = 0; ii < 3 * configs.J_coupled_max; ii = ii + 1)
        {
            phase_coupled(in, ii) = angle_trans * phase2[ii];
        }
    }

    bar.mark_as_completed();
    // Show cursor
    show_console_cursor(true);

    // write results in the output file.
    std::ofstream fp(configs.result_file());
    std::cout.precision(6);
    auto name_uncoupled = phases::partial_wave_name_uncoupled(configs);
    auto name_coupled = phases::partial_wave_name_coupled(configs);
    // Write partial wave names in the first line of the output file
    fp << std::fixed << "p_lab";
    for (const auto &name : name_uncoupled)
    {
        fp << "\t" << std::setw(15) << name;
    }
    for (const auto &name : name_coupled)
    {
        fp << "\t" << std::setw(15) << name;
    }
    fp << "\n";
    // Write phase shift values for each scattering angle and partial wave
    for (int ii = 0; ii < plablist.size(); ii = ii + 1)
    {
        fp << std::fixed << 1000.0 * plablist[ii];
        for (int jj = 0; jj < 2 * configs.J_uncoupled_max + 2; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_uncoupled(ii, jj);
        }
        for (int jj = 0; jj < 3 * configs.J_coupled_max; jj = jj + 1)
        {
            fp << "\t" << std::setw(15) << phase_coupled(ii, jj);
        }
        fp << '\n';
    }
    fp.close();
}

} // namespace evaluate

#endif