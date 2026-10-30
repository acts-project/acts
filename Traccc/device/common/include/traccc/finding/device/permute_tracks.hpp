// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/utils/logging.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <cassert>

namespace traccc::device {

/// Payload for the @c traccc::device::permute_tracks function
struct permute_tracks_payload {
  /**
   * @brief Number of track parameters (and links) to permute in this step
   */
  const unsigned int n_in_params;

  /**
   * @brief Index of the first link of the current step in the link vector
   */
  const unsigned int curr_links_idx;

  /**
   * @brief View object to the track parameters in their pre-permutation order
   */
  bound_track_parameters_collection_types::const_view in_params_view;

  /**
   * @brief View object to the parameter liveness values in their
   * pre-permutation order
   */
  vecmem::data::vector_view<const unsigned int> in_param_liveness_view;

  /**
   * @brief View object to the full vector of candidate links, of which the
   * current step's links are read in their pre-permutation order
   */
  vecmem::data::vector_view<const candidate_link> in_links_view;

  /**
   * @brief View object to the permutation to apply; element @c i holds the
   * pre-permutation index of the parameter that ends up at position @c i
   */
  vecmem::data::vector_view<const unsigned int> in_params_index_view;

  /**
   * @brief View object to the map from condensed parameter index to the slot
   * in the find_tracks scratch buffers (@c in_params_view and @c in_links_view)
   */
  vecmem::data::vector_view<const unsigned int> in_condensed_param_map_view;

  /**
   * @brief View object to the permuted track parameters (must not alias
   * @c in_params_view)
   */
  bound_track_parameters_collection_types::view out_params_view;

  /**
   * @brief View object to the permuted parameter liveness values (must not
   * alias @c in_param_liveness_view)
   */
  vecmem::data::vector_view<unsigned int> out_param_liveness_view;

  /**
   * @brief View object to a scratch vector receiving the current step's links
   * in permuted order; the caller copies it back into the link vector
   */
  vecmem::data::vector_view<candidate_link> out_links_view;

  /**
   * @brief View object to the per-link filtered parameters, written at the
   * permuted link index; may be empty if the MBF smoother is not used
   */
  bound_track_parameters_collection_types::view
      out_link_filtered_parameter_view;
};

/// Apply the sorted parameter order of a CKF step to the data itself
///
/// After sorting, the parameters, their liveness values and the current step's
/// links are gathered into sorted order, so that every subsequent kernel can
/// address them directly by thread index rather than through an index vector.
/// Because the links of this step move, this function also records the
/// filtered parameters per link at their final position.
///
/// @param[in] thread_id The index of the output slot to fill
/// @param[in] payload The payload for the function
///
TRACCC_HOST_DEVICE inline void permute_tracks(
    global_index_t thread_id, const permute_tracks_payload& payload) {
  if (thread_id >= payload.n_in_params) {
    return;
  }

  bound_track_parameters_collection_types::const_device in_params(
      payload.in_params_view);
  vecmem::device_vector<const unsigned int> in_param_liveness(
      payload.in_param_liveness_view);
  vecmem::device_vector<const candidate_link> in_links(payload.in_links_view);
  vecmem::device_vector<const unsigned int> in_params_index(
      payload.in_params_index_view);
  vecmem::device_vector<const unsigned int> in_condensed_param_map(
      payload.in_condensed_param_map_view);

  bound_track_parameters_collection_types::device out_params(
      payload.out_params_view);
  vecmem::device_vector<unsigned int> out_param_liveness(
      payload.out_param_liveness_view);
  vecmem::device_vector<candidate_link> out_links(payload.out_links_view);
  bound_track_parameters_collection_types::device out_link_filtered_parameter(
      payload.out_link_filtered_parameter_view);

  const unsigned int param_id = in_params_index.at(thread_id);

  const unsigned int uncondensed_id = in_condensed_param_map.at(param_id);

  out_params.at(thread_id) = in_params.at(uncondensed_id);
  out_param_liveness.at(thread_id) = in_param_liveness.at(param_id);
  out_links.at(thread_id) = in_links.at(uncondensed_id);
  if (out_link_filtered_parameter.capacity() > 0) {
    out_link_filtered_parameter.at(payload.curr_links_idx + thread_id) =
        in_params.at(uncondensed_id);
  }
}

}  // namespace traccc::device
