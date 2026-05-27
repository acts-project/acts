/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/csv/cell.hpp"
#include "traccc/io/csv/hit.hpp"
#include "traccc/io/csv/measurement.hpp"
#include "traccc/io/csv/measurement_hit_id.hpp"
#include "traccc/io/csv/particle.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/utils/seed_generator.hpp"

// Vecmem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <map>
#include <string>

namespace traccc {

struct event_data {

    public:
    event_data() = delete;

    /// Event data constructor
    ///
    /// @param[in] event_dir Event data directory
    /// data
    /// @param[in] event_id  Event id
    /// @param[in] use_acts_geom_source  Use acts geometry source
    /// @param[in] det       detray detector
    /// @param[in] format    file format
    /// @param[in] include_silicon_cells Use silicon cell data in object
    /// construction
    ///
    event_data(const std::string& event_dir, const std::size_t event_id,
               vecmem::memory_resource& resource,
               bool use_acts_geom_source = false,
               const host_detector* det = nullptr,
               data_format format = data_format::csv,
               bool include_silicon_cells = false);

    /// Fill the member variables related to CCA
    ///
    /// @param[in] cells cell EDM
    /// @param[in] cca_clusters cluster EDM from CCL algorithm
    /// @param[in] cca_measurements measurement EDM from measurement creation
    /// @param[in] det_desc    Detector design description
    ///
    void fill_cca_result(
        const edm::silicon_cell_collection::host& cells,
        const edm::silicon_cluster_collection::host& cca_clusters,
        const edm::measurement_collection::host& cca_measurements,
        const detector_conditions_description::host& det_cond);

    /// Generate truth candidate used for truth fitting
    ///
    /// @param[out] truth_candidates Truth candidates
    /// @param[in] sg Seed generator for fitting
    /// @param[in] resource vecmem memory resource
    ///
    template <typename detector_type>
    void generate_truth_candidates(
        edm::track_container<default_algebra>::host& truth_candidates,
        edm::measurement_collection::host& truth_measurements,
        seed_generator<detector_type>& sg, vecmem::memory_resource& resource,
        float pt_cut = 0.f) {
        for (auto const& [ptc, measurements] : m_ptc_to_meas_map) {

            const auto& param = m_meas_to_param_map.at(measurements[0]);
            const free_track_parameters<> free_param(param.first, 0.f,
                                                     param.second, ptc.charge);

            auto ptc_particle =
                detail::particle_from_pdg_number<scalar>(ptc.particle_type);

            if (ptc_particle.pdg_num() == 0) {
                // TODO: Add some debug logging here.
                continue;
            } else if (free_param.pT(ptc_particle.charge()) <= pt_cut) {
                continue;
            }

            auto seed_params =
                sg(measurements[0].surface_link(), free_param, ptc_particle);

            // Record the measurements, and remember their indices.
            vecmem::vector<edm::track_constituent_link> meas_links{&resource};
            truth_measurements.reserve(truth_measurements.size() +
                                       measurements.size());
            meas_links.reserve(measurements.size());
            for (const auto& meas : measurements) {
                meas_links.push_back(
                    {edm::track_constituent_link::measurement,
                     static_cast<unsigned int>(truth_measurements.size())});
                truth_measurements.push_back(meas);
            }

            // Record the truth track candidate.
            edm::track_collection<traccc::default_algebra>::host::object_type
                track;
            track.params() = seed_params;
            track.constituent_links() = meas_links;
            truth_candidates.tracks.push_back(track);
        }
    }

    /// All internally defined measurements
    edm::measurement_collection::host m_measurements;
    /// Proxy type for the following maps
    using measurement_proxy = edm::measurement_collection::host::object_type;

    // Measurement map
    std::map<measurement_id_type, measurement_proxy> m_measurement_map;
    // Particle map
    std::map<particle_id, particle> m_particle_map;
    // Measurement to the contributing particle map
    std::map<measurement_proxy, std::map<particle, std::size_t>>
        m_meas_to_ptc_map;
    // CCA measurement to the contributing particle map
    std::map<measurement_proxy, std::map<particle, std::size_t>>
        m_found_meas_to_ptc_map;
    // Particle to its Measurements map
    std::map<particle, std::vector<measurement_proxy>> m_ptc_to_meas_map;
    // Measurement to its track parameter map
    std::map<measurement_proxy, std::pair<point3, point3>> m_meas_to_param_map;
    // CCA measurement to its track parameter map
    std::map<measurement_proxy, std::pair<point3, point3>>
        m_found_meas_to_param_map;
    // Cell to particle map
    std::map<io::csv::cell, particle> m_cell_to_particle_map;

    // Input arguments
    const std::string m_event_dir;
    const std::size_t m_event_id;
    std::reference_wrapper<vecmem::memory_resource> m_mr;

    private:
    void setup_csv(bool use_acts_geom_source, const host_detector* det,
                   bool include_silicon_cells);
};

}  // namespace traccc
