/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/options/track_matching.hpp"

#include <format>

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"

namespace traccc::opts {

track_matching::track_matching() : interface("Track Matching Options") {
    m_desc.add_options()("track-matching-ratio",
                         boost::program_options::value(&m_matching_ratio)
                             ->default_value(m_matching_ratio),
                         "Minimum track state matching ratio");
    m_desc.add_options()("track-double-matching",
                         boost::program_options::value(&m_double_matching)
                             ->default_value(m_double_matching),
                         "Enable double truth-reco matching");
}

track_matching::operator track_matching_config() const {
    return track_matching_config{.matching_ratio = m_matching_ratio,
                                 .double_matching = m_double_matching};
}

std::unique_ptr<configuration_printable> track_matching::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Matching ratio", std::format("{}", m_matching_ratio)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Double matching", std::format("{}", m_double_matching)));

    return cat;
}
}  // namespace traccc::opts
