/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/fitting/triplet_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"
#include "traccc/utils/propagation.hpp"
// GBTS include for dummy input (not implemented)
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>
#include <memory>

namespace traccc {

/// Algorithm performing the full chain of track reconstruction
///
/// At least as much as is implemented in the project at any given moment.
///
class full_chain_algorithm
    : public algorithm<edm::track_collection<default_algebra>::host(
          const edm::silicon_cell_collection::host&)>,
      public messaging {

    public:
    /// @name Type declaration(s)
    /// @{

    /// Clusterization algorithm type
    using clustering_algorithm = host::clusterization_algorithm;
    /// Spacepoint formation algorithm type
    using spacepoint_formation_algorithm =
        traccc::host::silicon_pixel_spacepoint_formation_algorithm;
    /// Track finding algorithm type
    using finding_algorithm =
        traccc::host::combinatorial_kalman_filter_algorithm;
    /// Track fitting algorithm type
    using fitting_algorithm = traccc::host::kalman_fitting_algorithm;

    /// @}

    /// Algorithm constructor
    ///
    /// @param mr The memory resource to use for the intermediate and result
    ///           objects
    /// @param dummy This is not used anywhere. Allows templating CPU/Device
    /// algorithm.
    ///
    full_chain_algorithm(
        vecmem::memory_resource& mr,
        const clustering_algorithm::config_type& dummy,
        const seedfinder_config& finder_config,
        const spacepoint_grid_config& grid_config,
        const seedfilter_config& filter_config,
        const gbts_seedfinder_config& gbts_config,
        const track_params_estimation_config& track_params_estimation_config,
        const finding_algorithm::config_type& finding_config,
        const fitting_algorithm::config_type& fitting_config,
        const detector_design_description::host& det_descr,
        const detector_conditions_description::host& det_cond,
        const magnetic_field& field, const host_detector* detector,
        std::unique_ptr<const traccc::Logger> logger,
        const bool useGBTS = false);

    /// Reconstruct track parameters in the entire detector
    ///
    /// @param cells The cells for every detector module in the event
    /// @return The track parameters reconstructed
    ///
    output_type operator()(
        const edm::silicon_cell_collection::host& cells) const override;

    /// Reconstruct track seeds in the entire detector
    ///
    /// @param cells The cells for every detector module in the event
    /// @return The track seeds reconstructed
    ///
    bound_track_parameters_collection_types::host seeding(
        const edm::silicon_cell_collection::host& cells) const;

    private:
    /// Memory resource
    std::reference_wrapper<vecmem::memory_resource> m_mr;
    /// Vecmem copy object
    std::unique_ptr<vecmem::copy> m_copy;
    /// Constant B field for the (seed) track parameter estimation
    traccc::vector3 m_field_vec;
    /// Constant B field for the track finding and fitting
    magnetic_field m_field;

    /// Detector description
    std::reference_wrapper<const detector_design_description::host> m_det_descr;
    std::reference_wrapper<const detector_conditions_description::host>
        m_det_cond;
    /// Detector
    const host_detector* m_detector;

    /// @name Sub-algorithms used by this full-chain algorithm
    /// @{

    /// Clusterization algorithm
    clustering_algorithm m_clusterization;
    /// Spacepoint formation algorithm
    spacepoint_formation_algorithm m_spacepoint_formation;
    /// Seeding algorithm
    host::seeding_algorithm m_seeding;
    /// Track parameter estimation algorithm
    host::track_params_estimation m_track_parameter_estimation;

    /// Track finding algorithm
    finding_algorithm m_finding;
    /// Track fitting algorithm
    fitting_algorithm m_fitting;

    /// @}

    /// @name Algorithm configurations
    /// @{

    /// Configuration for the seed finding
    seedfinder_config m_finder_config;
    /// Configuration for the spacepoint grid formation
    spacepoint_grid_config m_grid_config;
    /// Configuration for the seed filtering
    seedfilter_config m_filter_config;
    /// placeholder GBTS config
    [[maybe_unused]] gbts_seedfinder_config m_gbts_config;
    /// Configuration for track parameter estimation
    track_params_estimation_config m_track_params_estimation_config;

    /// Configuration for the track finding
    finding_algorithm::config_type m_finding_config;
    /// Configuration for the track fitting
    fitting_algorithm::config_type m_fitting_config;

    const bool usingGBTS;

    /// @}
};  // class full_chain_algorithm

}  // namespace traccc
