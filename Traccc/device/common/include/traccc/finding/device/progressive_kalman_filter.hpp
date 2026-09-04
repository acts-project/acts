/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/actors/measurement_kalman_updater.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/track_state_candidate.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::progressive_kalman_filter
/// function
struct progressive_kalman_filter_payload {
  using algebra_t = default_algebra;
  using scalar_t = detray::dscalar<algebra_t>;

  /**
   * @brief View object to the vector of bound track parameters
   */
  bound_track_parameters_collection_types::view seeds_view;

  /**
   * @brief View object to the vector of input measurements
   *
   * @warning Measurements on the same surface must be adjacent
   */
  edm::measurement_collection::const_view measurements_view;

  /**
   * @brief View object to the vector of measurement index ranges per surface
   */
  vecmem::data::vector_view<unsigned int> measurement_ranges_view;

  /**
   * @brief View object to the collected track statistic data
   */
  vecmem::data::vector_view<track_stats<scalar_t>> track_stats_view;

  /**
   * @brief View object to the collected track state data
   */
  vecmem::data::vector_view<track_state_candidate> track_cand_view;
  vecmem::data::vector_view<filtered_track_state_candidate<algebra_t>>
      filtered_track_cand_view;
  vecmem::data::vector_view<full_track_state_candidate<algebra_t>>
      full_track_cand_view;

  /**
   * @brief View object to the vector of track candidates
   */
  edm::track_container<algebra_t>::view tracks_view;
};

/// Function that runs Kalman filter based track following
///
/// @param[in] globalIndex        The index of the current thread
/// @param[in] cfg                Track finding config object
/// @param[inout] payload         The function call payload
///
template <typename propagator_t>
TRACCC_HOST_DEVICE inline void progressive_kalman_filter(
    global_index_t globalIndex, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const typename propagator_t::stepper_type::magnetic_field_type& field_data,
    vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>
        surfaces_view,
    const progressive_kalman_filter_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/progressive_kalman_filter.ipp"
