/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/telescope_detector.hpp"

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <sstream>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

telescope_detector::telescope_detector()
    : interface("Telescope Detector Options") {

    m_desc.add_options()("empty-material", po::bool_switch(&empty_material),
                         "Build detector without materials");
    m_desc.add_options()("n-planes",
                         po::value(&n_planes)->default_value(n_planes),
                         "Number of planes");
    m_desc.add_options()("thickness-mm",
                         po::value(&thickness)->default_value(thickness),
                         "Slab thickness in [mm]");
    m_desc.add_options()("spacing", po::value(&spacing)->default_value(spacing),
                         "Space between planes in [mm]");
    m_desc.add_options()("smearing-um",
                         po::value(&smearing)->default_value(smearing),
                         "Measurement smearing in [um]");
    m_desc.add_options()("half-length-mm",
                         po::value(&half_length)->default_value(half_length),
                         "Half length of plane [mm]");
    m_desc.add_options()("align-vector",
                         po::value(&align_vector)
                             ->value_name("X:Y:Z")
                             ->default_value(align_vector),
                         "Vector for plane placement");
}

void telescope_detector::read(const po::variables_map &) {

    thickness *= traccc::unit<float>::mm;
    spacing *= traccc::unit<float>::mm;
    smearing *= traccc::unit<float>::um;
    half_length *= traccc::unit<float>::mm;
}

std::unique_ptr<configuration_printable> telescope_detector::as_printable()
    const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Empty material", empty_material ? "yes" : "no"));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of planes", std::to_string(n_planes)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Slab thickness",
        std::to_string(thickness / traccc::unit<float>::mm) + " mm"));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Spacing", std::to_string(spacing / traccc::unit<float>::mm) + " mm"));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Smearing",
        std::to_string(thickness / traccc::unit<float>::um) + " um"));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Half length",
        std::to_string(half_length / traccc::unit<float>::mm) + " mm"));
    std::ostringstream align_ss;
    align_ss << align_vector;
    cat->add_child(std::make_unique<configuration_kv_pair>("Alignment axis",
                                                           align_ss.str()));

    return cat;
}

}  // namespace traccc::opts
