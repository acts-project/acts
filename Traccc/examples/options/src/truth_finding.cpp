/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/options/truth_finding.hpp"

#include <format>

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"

namespace traccc::opts {

truth_finding::truth_finding() : interface("Truth Track Finding Options") {
    m_desc.add_options()(
        "truth-finding-min-pt",
        boost::program_options::value(&m_pT_min)->default_value(m_pT_min),
        "Candidate particle pT cut [GeV]");
    m_desc.add_options()(
        "truth-finding-min-z",
        boost::program_options::value(&m_z_min)->default_value(m_z_min),
        "Candidate particle min z cut [mm]");
    m_desc.add_options()(
        "truth-finding-max-z",
        boost::program_options::value(&m_z_max)->default_value(m_z_max),
        "Candidate particle max z cut [mm]");
    m_desc.add_options()(
        "truth-finding-max-r",
        boost::program_options::value(&m_r_max)->default_value(m_r_max),
        "Candidate particle max r cut [mm]");
    m_desc.add_options()(
        "truth-finding-max-eta",
        boost::program_options::value(&m_eta_max)->default_value(m_eta_max),
        "Candidate particle max eta cut");
    m_desc.add_options()("truth-finding-process-id",
                         boost::program_options::value(&m_process_id)
                             ->default_value(m_process_id),
                         "Candidate particle is from a selected process");
    m_desc.add_options()("truth-finding-min-track-candidates",
                         boost::program_options::value(&m_min_track_candidates)
                             ->default_value(m_min_track_candidates),
                         "Minimum track candidates on track");
}

void truth_finding::read(const boost::program_options::variables_map &) {
    m_pT_min *= unit<float>::GeV;
    m_z_min *= unit<float>::mm;
    m_z_max *= unit<float>::mm;
    m_r_max *= unit<float>::mm;
}

truth_finding::operator truth_matching_config() const {
    return truth_matching_config{
        .pT_min = m_pT_min,
        .z_min = m_z_min,
        .z_max = m_z_max,
        .r_max = m_r_max,
        .eta_max = m_eta_max,
        .process_id = m_process_id,
        .min_track_candidates = m_min_track_candidates};
}

std::unique_ptr<configuration_printable> truth_finding::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum pT", std::format("{} GeV", m_pT_min / unit<float>::GeV)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum z", std::format("{} mm", m_z_min / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum z", std::format("{} mm", m_z_max / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum r", std::format("{} mm", m_r_max / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum eta", std::format("{}", m_eta_max)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Process ID", std::format("{}", m_process_id)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum track candidates", std::format("{}", m_min_track_candidates)));

    return cat;
}
}  // namespace traccc::opts
