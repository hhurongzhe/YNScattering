add_rules("mode.release")
add_requires("openmp")

target("YNScattering")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("aggressive")
    add_packages("openmp")
    add_files("main.cpp")


