# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

traccc_add_executable(convert_csv_bfield "convert_csv_bfield.cpp"
    LINK_LIBRARIES traccc::core covfie::core Boost::program_options)

traccc_add_executable(generate_constant_bfield "generate_constant_bfield.cpp"
    LINK_LIBRARIES traccc::core covfie::core Boost::program_options)
