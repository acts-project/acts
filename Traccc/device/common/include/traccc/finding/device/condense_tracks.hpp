/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

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

/// (Event Data) Payload for the @c traccc::device::condense_tracks function
struct condense_tracks_payload {
  /**
   * @brief The total number of input parameters
   */
  const unsigned int n_in_params;
  const unsigned int step;
  const unsigned int curr_links_idx;
  const unsigned int max_num_branches_per_surface;
  const unsigned int min_track_candidates_per_track;
  const unsigned int max_track_candidates_per_track;

  /**
   * @brief View object to the link vector
   */
  vecmem::data::vector_view<candidate_link> links_view;

  /**
   * @brief View object to the vector of track parameters
   */
  bound_track_parameters_collection_types::const_view in_params_view;

  /**
   * @brief View object to the temporary track parameter vector
   */
  bound_track_parameters_collection_types::const_view in_tmp_params_view;

  /**
   * @brief View object to the temporary link vector
   */
  vecmem::data::vector_view<const candidate_link> in_tmp_links_view;

  vecmem::data::vector_view<const unsigned int> in_params_index_view;

  /**
   * @brief View object to the output track parameter vector
   */
  bound_track_parameters_collection_types::view out_params_view;

  /**
   * @brief View object to the output track parameter liveness vector
   */
  vecmem::data::vector_view<unsigned int> out_params_liveness_view;

  /**
   * @brief View object to the vector of tips
   */
  vecmem::data::vector_view<unsigned int> tips_view;

  /**
   * @brief Vector to hold the number of track states per tip
   */
  vecmem::data::vector_view<unsigned int> tip_lengths_view;

  vecmem::data::vector_view<bound_matrix<default_algebra>> jacobian_view;
  vecmem::data::vector_view<bound_matrix<default_algebra>> tmp_jacobian_view;
  bound_track_parameters_collection_types::view link_predicted_parameter_view;
  bound_track_parameters_collection_types::view link_filtered_parameter_view;
};

/// Condense the per-input-parameter outputs of @c traccc::device::find_tracks
/// into a tightly-packed array using the prefix-summed
/// @c in_params_index_view, then write the corresponding Jacobians, the
/// optional predicted/filtered parameters, and emit tips on the last step.
///
/// @param[in] thread_id The index of the current thread
/// @param[inout] payload The function call payload
///
TRACCC_HOST_DEVICE inline void condense_tracks(
    global_index_t thread_id, const condense_tracks_payload& payload) {
  const unsigned int in_param_id = thread_id;

  if (in_param_id >= payload.n_in_params) {
    return;
  }

  const bool last_step =
      payload.step == payload.max_track_candidates_per_track - 1;

  vecmem::device_vector<const unsigned int> in_params_index(
      payload.in_params_index_view);
  vecmem::device_vector<const candidate_link> in_tmp_links(
      payload.in_tmp_links_view);
  vecmem::device_vector<candidate_link> links(payload.links_view);
  bound_track_parameters_collection_types::const_device in_params(
      payload.in_params_view);
  bound_track_parameters_collection_types::const_device in_tmp_params(
      payload.in_tmp_params_view);
  bound_track_parameters_collection_types::device out_params(
      payload.out_params_view);
  vecmem::device_vector<unsigned int> out_params_liveness(
      payload.out_params_liveness_view);
  vecmem::device_vector<unsigned int> tips(payload.tips_view);
  vecmem::device_vector<unsigned int> tip_lengths(payload.tip_lengths_view);
  vecmem::device_vector<bound_matrix<default_algebra>> jacobian(
      payload.jacobian_view);
  vecmem::device_vector<bound_matrix<default_algebra>> tmp_jacobian(
      payload.tmp_jacobian_view);
  bound_track_parameters_collection_types::device link_predicted_parameters(
      payload.link_predicted_parameter_view);
  bound_track_parameters_collection_types::device link_filtered_parameters(
      payload.link_filtered_parameter_view);

  unsigned int num_parameters;
  unsigned int idx;

  if (in_param_id == 0) {
    idx = 0;
    num_parameters = in_params_index.at(in_param_id);
  } else {
    idx = in_params_index.at(in_param_id - 1);
    num_parameters = in_params_index.at(in_param_id) - idx;
  }

  for (unsigned int i = 0; i < num_parameters; ++i) {
    const unsigned int in_offset =
        in_param_id * payload.max_num_branches_per_surface + i;
    const unsigned int param_out_index = idx + i;
    const unsigned int link_out_index =
        param_out_index + payload.curr_links_idx;

    out_params.at(param_out_index) = in_tmp_params.at(in_offset);
    out_params_liveness.at(param_out_index) =
        static_cast<unsigned int>(!last_step);
    links.at(link_out_index) = in_tmp_links.at(in_offset);

    if (payload.tmp_jacobian_view.ptr() != nullptr) {
      assert(payload.jacobian_view.ptr() != nullptr);
      jacobian.at(link_out_index) = tmp_jacobian.at(in_param_id);
    }

    if (link_filtered_parameters.capacity() > 0) {
      link_filtered_parameters.at(link_out_index) = in_tmp_params.at(in_offset);
    }

    if (link_predicted_parameters.capacity() > 0) {
      link_predicted_parameters.at(link_out_index) = in_params.at(in_param_id);
    }

    const unsigned int n_cands =
        payload.step + 1 - in_tmp_links.at(in_offset).n_skipped;

    if (last_step && n_cands >= payload.min_track_candidates_per_track) {
      TRACCC_ERROR_DEVICE("Create tip: Max no. candidates");
      auto tip_pos = tips.push_back(link_out_index);
      tip_lengths.at(tip_pos) = n_cands;
    }
  }
}

}  // namespace traccc::device
