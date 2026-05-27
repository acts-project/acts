/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/detector.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

detector::detector() : interface("Detector Options") {

    namespace po = boost::program_options;

    m_desc.add_options()(
        "detector-file",
        po::value(&detector_file)->default_value(detector_file),
        "Detector file");
    m_desc.add_options()(
        "material-file",
        po::value(&material_file)->default_value(material_file),
        "Material file");
    m_desc.add_options()("grid-file",
                         po::value(&grid_file)->default_value(grid_file),
                         "Surface grid file");
    m_desc.add_options()(
        "digitization-file",
        po::value(&digitization_file)->default_value(digitization_file),
        "Digitization file");
    m_desc.add_options()(
        "conditions-file",
        po::value(&conditions_file)->default_value(conditions_file),
        "Conditions file");
}

std::unique_ptr<configuration_printable> detector::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>("Detector file",
                                                           detector_file));
    cat->add_child(std::make_unique<configuration_kv_pair>("Material file",
                                                           material_file));
    cat->add_child(std::make_unique<configuration_kv_pair>("Surface grid file",
                                                           grid_file));
    cat->add_child(std::make_unique<configuration_kv_pair>("Digitization file",
                                                           digitization_file));
    cat->add_child(std::make_unique<configuration_kv_pair>("Conditions file",
                                                           conditions_file));

    return cat;
}

}  // namespace traccc::opts
