/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/options/seed_matching.hpp"

#include <format>

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"

namespace traccc::opts {

seed_matching::seed_matching() : interface("Seed Matching Options") {
    m_desc.add_options()("seed-matching-ratio",
                         boost::program_options::value(&m_matching_ratio)
                             ->default_value(m_matching_ratio),
                         "Minimum track state matching ratio");
}

seed_matching::operator seed_matching_config() const {
    return seed_matching_config{
        .matching_ratio = m_matching_ratio,
    };
}

std::unique_ptr<configuration_printable> seed_matching::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Matching ratio", std::format("{}", m_matching_ratio)));

    return cat;
}
}  // namespace traccc::opts
