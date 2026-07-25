/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::build_tracks function
struct build_tracks_payload {
  /**
   * @brief View object to the vector of measurements
   */
  bound_track_parameters_collection_types::const_view seeds_view;

  /**
   * @brief View object to the vector of candidate links
   */
  vecmem::data::vector_view<const candidate_link> links_view;

  /**
   * @brief View object to the vector of tips
   */
  vecmem::data::vector_view<const unsigned int> tips_view;

  /**
   * @brief View object to the vector of track candidates
   */
  edm::track_container<default_algebra>::view tracks_view;

  /**
   * @brief Optional mapping from tip index to output index
   */
  vecmem::data::vector_view<const unsigned int> tip_to_output_map;
  bound_matrix<default_algebra>* jacobian_ptr = nullptr;
  bound_track_parameters_collection_types::view link_predicted_parameter_view;
  bound_track_parameters_collection_types::view link_filtered_parameter_view;
};

/// Function for building full tracks from the link container:
/// The full tracks are built using the link container and tip link container.
/// Since every link holds an information of the link from the previous step,
/// we can build a full track by iterating from a tip link backwardly.
///
/// @param[in] globalIndex         The index of the current thread
/// @param[in] cfg                    Track finding config object
/// @param[inout] payload      The function call payload
///
TRACCC_HOST_DEVICE inline void build_tracks(
    global_index_t globalIndex, bool run_mbf,
    const measurement_selector::config calib_cfg,
    const build_tracks_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/build_tracks.ipp"
