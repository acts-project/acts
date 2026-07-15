/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/efficiency/seeding_performance_writer.hpp"

#include "duplication_plot_tool.hpp"
#include "eff_plot_tool.hpp"
#include "track_classification.hpp"

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <TFile.h>
#endif  // TRACCC_HAVE_ROOT

// System include(s).
#include <iostream>
#include <memory>
#include <stdexcept>

namespace traccc {
namespace details {

struct seeding_performance_writer_data {

    /// Constructor
    seeding_performance_writer_data(
        const seeding_performance_writer::config& cfg)
        : m_eff_plot_tool({cfg.var_binning}),
          m_duplication_plot_tool({cfg.var_binning}) {}

    /// Plot tool for efficiency
    eff_plot_tool m_eff_plot_tool;
    eff_plot_tool::eff_plot_cache m_eff_plot_cache;

    /// Plot tool for duplication rate
    duplication_plot_tool m_duplication_plot_tool;
    duplication_plot_tool::duplication_plot_cache m_duplication_plot_cache;

    std::map<event_data::measurement_proxy, std::map<particle, std::size_t>>
        m_measurement_particle_map;
    std::map<std::uint64_t, particle> m_particle_map;

};  // struct seeding_performance_writer_data

}  // namespace details

seeding_performance_writer::seeding_performance_writer(
    const config& cfg, std::unique_ptr<const traccc::Logger> logger)
    : messaging(std::move(logger)),
      m_cfg(cfg),
      m_data(std::make_unique<details::seeding_performance_writer_data>(cfg)) {

    m_data->m_eff_plot_tool.book("seeding", m_data->m_eff_plot_cache);
    m_data->m_duplication_plot_tool.book("seeding",
                                         m_data->m_duplication_plot_cache);
}

seeding_performance_writer::~seeding_performance_writer() {}

void seeding_performance_writer::write(
    const edm::seed_collection::const_view& seeds_view,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const edm::measurement_collection::const_view& measurements_view,
    const event_data& evt_data) {

    std::map<particle_id, std::size_t> match_counter;

    // Create the device collections.
    const edm::seed_collection::const_device seeds(seeds_view);
    const edm::spacepoint_collection::const_device spacepoints(
        spacepoints_view);
    const edm::measurement_collection::const_device measurements(
        measurements_view);

    const std::size_t n_seeds = seeds.size();

    std::size_t total_fake_seeds = 0;

    // Iterate over the seeds.
    for (edm::seed_collection::const_device::size_type i = 0u; i < seeds.size();
         ++i) {
        const auto sd = seeds.at(i);

        std::vector<particle_hit_count> particle_hit_counts;

        if (sd.bottom_index() >= spacepoints.size() ||
            sd.middle_index() >= spacepoints.size() ||
            sd.top_index() >= spacepoints.size()) {
            TRACCC_WARNING("Seed with index "
                           << i << " (out of " << seeds.size()
                           << ") has invalid spacepoint (B: "
                           << sd.bottom_index() << ", M: " << sd.middle_index()
                           << ", T: " << sd.top_index()
                           << ", #SP: " << spacepoints.size() << ")");
            continue;
        }

        if (spacepoints.at(sd.bottom_index()).measurement_index_1() >=
                measurements.size() ||
            spacepoints.at(sd.middle_index()).measurement_index_1() >=
                measurements.size() ||
            spacepoints.at(sd.top_index()).measurement_index_1() >=
                measurements.size()) {
            TRACCC_WARNING(
                "Seed with index "
                << i << " (out of " << seeds.size()
                << ") has invalid measurement (B: "
                << spacepoints.at(sd.bottom_index()).measurement_index_1()
                << ", M: "
                << spacepoints.at(sd.middle_index()).measurement_index_1()
                << ", T: "
                << spacepoints.at(sd.top_index()).measurement_index_1()
                << ", #SP: " << spacepoints.size() << ")");
            continue;
        }

        // Get the measurements for this seed.
        assert(spacepoints.at(sd.bottom_index()).measurement_index_2() ==
               edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX);
        assert(spacepoints.at(sd.middle_index()).measurement_index_2() ==
               edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX);
        assert(spacepoints.at(sd.top_index()).measurement_index_2() ==
               edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX);
        std::array<event_data::measurement_proxy, 3> seed_measurements{
            measurements.at(
                spacepoints.at(sd.bottom_index()).measurement_index_1()),
            measurements.at(
                spacepoints.at(sd.middle_index()).measurement_index_1()),
            measurements.at(
                spacepoints.at(sd.top_index()).measurement_index_1())};

        if (!evt_data.m_found_meas_to_ptc_map.empty()) {
            particle_hit_counts = identify_contributing_particles(
                seed_measurements, evt_data.m_found_meas_to_ptc_map);
        } else {
            particle_hit_counts = identify_contributing_particles(
                seed_measurements, evt_data.m_meas_to_ptc_map);
        }

        // Consider it being matched if hit counts is larger than the half
        // of the number of measurements
        assert(seed_measurements.size() > 0u);
        if (static_cast<double>(particle_hit_counts.at(0).hit_counts) /
                static_cast<double>(seed_measurements.size()) >
            m_cfg.seed_truth_config.matching_ratio) {
            auto pid = particle_hit_counts.at(0).ptc.particle_id;
            match_counter[pid]++;
        } else {
            total_fake_seeds++;
        }
    }

    std::size_t total_truth_particles = 0;
    std::size_t total_matched_truth_particles = 0;
    std::size_t total_duplicate_seeds = 0;

    for (auto const& [pid, ptc] : evt_data.m_particle_map) {
        std::size_t num_measurements = 0;

        // Find the number of measurements belonging to this track
        if (auto it = evt_data.m_ptc_to_meas_map.find(ptc);
            it != evt_data.m_ptc_to_meas_map.end()) {
            num_measurements = it->second.size();
        } else {
            continue;
        }

        // Count only charged particles which satisfiy pT_cut
        if (ptc.charge == 0 ||
            vector::perp(ptc.momentum) < m_cfg.truth_config.pT_min ||
            ptc.vertex[2] < m_cfg.truth_config.z_min ||
            ptc.vertex[2] > m_cfg.truth_config.z_max ||
            vector::perp(ptc.vertex) > m_cfg.truth_config.r_max ||
            std::abs(vector::eta(ptc.momentum)) > m_cfg.truth_config.eta_max ||
            (m_cfg.truth_config.process_id >= 0 &&
             m_cfg.truth_config.process_id != ptc.process) ||
            num_measurements < m_cfg.truth_config.min_track_candidates) {
            continue;
        }

        total_truth_particles++;

        bool is_matched = false;
        std::size_t n_matched_seeds_for_particle = 0;
        auto it = match_counter.find(pid);
        if (it != match_counter.end()) {
            is_matched = true;
            n_matched_seeds_for_particle = it->second;
            total_matched_truth_particles++;
            total_duplicate_seeds += n_matched_seeds_for_particle - 1;
        }

        m_data->m_eff_plot_tool.fill(m_data->m_eff_plot_cache, ptc, is_matched);
        m_data->m_duplication_plot_tool.fill(m_data->m_duplication_plot_cache,
                                             ptc,
                                             n_matched_seeds_for_particle - 1);
    }

    TRACCC_INFO("Total number of truth particles was "
                << total_truth_particles);
    TRACCC_INFO("Total number of found seeds was " << n_seeds);
    TRACCC_INFO("Total number of seed-matched particles was "
                << total_matched_truth_particles);
    TRACCC_INFO("Total number of duplicate seeds was "
                << total_duplicate_seeds);
    TRACCC_INFO("Total number of fake seeds was " << total_fake_seeds);
    TRACCC_INFO("Total seed efficiency was "
                << (100. * static_cast<double>(total_matched_truth_particles) /
                    static_cast<double>(total_truth_particles))
                << "%");
    TRACCC_INFO("Total seed duplicate rate was "
                << (static_cast<double>(total_duplicate_seeds) /
                    static_cast<double>(total_matched_truth_particles)));
    TRACCC_INFO("Total seed fake rate was "
                << (static_cast<double>(total_fake_seeds) /
                    static_cast<double>(total_truth_particles)));
}

void seeding_performance_writer::finalize() {

#ifdef TRACCC_HAVE_ROOT
    // Open the output file.
    std::unique_ptr<TFile> ofile(
        TFile::Open(m_cfg.file_path.c_str(), m_cfg.file_mode.c_str()));
    if ((!ofile) || ofile->IsZombie()) {
        throw std::runtime_error("Could not open output file \"" +
                                 m_cfg.file_path + "\" in mode \"" +
                                 m_cfg.file_mode + "\"");
    }
    ofile->cd();
#else
    std::cout << "ROOT file \"" << m_cfg.file_path << "\" is NOT created"
              << std::endl;
#endif  // TRACCC_HAVE_ROOT

    m_data->m_eff_plot_tool.write(m_data->m_eff_plot_cache);
    m_data->m_duplication_plot_tool.write(m_data->m_duplication_plot_cache);
}

}  // namespace traccc
