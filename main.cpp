#include "include/evaluate.hpp"
#include "include/lib_define.hpp"
#include "include/observable.hpp"

int main()
{
    // config file initializing.
    auto ini = inifile_system::inifile("YN.ini");
    if (!ini.good())
    {
        std::cerr << ini.error() << std::endl;
        exit(-1);
    }
    auto configs = YN::YN_configs(ini);
    std::cout << configs.result_file() << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    //* main program:
    // evaluate::evaluate_phase_unmatched(configs);

    evaluate::evaluate_phase_matched(configs);

    // observable::differential_cross_section(configs);

    // observable::total_cross_section(configs);

    // observable::differential_cross_section_multichannel(configs);

    // observable::total_cross_section_multichannel(configs);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout.precision(4);
    std::cout << "Duration: " << duration.count() << " seconds\n";
}