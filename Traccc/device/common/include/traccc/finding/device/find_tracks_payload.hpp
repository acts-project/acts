/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <cstdint>
#include <utility>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::find_tracks function
struct find_tracks_payload {
  /**
   * @brief View object to the vector of bound track parameters
   *
   * @warning Measurements on the same surface must be adjacent
   */
  edm::measurement_collection::const_view measurements_view;

  /**
   * @brief View object to the vector of track parameters
   */
  bound_track_parameters_collection_types::const_view in_params_view;

  /**
   * @brief View object to the vector of boolean-like integers describing the
   * liveness of each parameter
   */
  vecmem::data::vector_view<const unsigned int> in_params_liveness_view;

  /**
   * @brief The total number of input parameters
   */
  unsigned int n_in_params;

  /**
   * @brief View object to the vector of measurement index ranges per surface
   */
  vecmem::data::vector_view<const unsigned int> measurement_ranges_view;

  /**
   * @brief View object to the link vector
   */
  vecmem::data::vector_view<candidate_link> links_view;

  /**
   * @brief Index in the link vector at which the previous step starts
   */
  const unsigned int prev_links_idx;

  /**
   * @brief The current step identifier
   */
  unsigned int step;

  /**
   * @brief View object to the output parameter counting vector
   */
  vecmem::data::vector_view<unsigned int> out_params_per_in_param_view;

  /**
   * @brief View object to the vector of tips
   */
  vecmem::data::vector_view<unsigned int> tips_view;

  /**
   * @brief Vector to hold the number of track states per tip
   */
  vecmem::data::vector_view<unsigned int> tip_lengths_view;

  /**
   * @brief View object to the vector of the number of tracks per initial
   * input seed
   */
  vecmem::data::vector_view<unsigned int> n_tracks_per_seed_view;

  /**
   * @brief View object to the temporary track parameter vector
   */
  bound_track_parameters_collection_types::view tmp_params_view;

  /**
   * @brief View object to the temporary link vector
   */
  vecmem::data::vector_view<candidate_link> tmp_links_view;
};

}  // namespace traccc::device
