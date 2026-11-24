// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void select_hits(const scalar* values, const scalar mode,
                        const unsigned int n_hits,
                        const vecmem::device_vector<const global_index_t> hits,
                        vecmem::device_vector<global_index_t> filtered_hits) {
  static constexpr scalar tolerance = 0.001f;
  for (unsigned int i = 0; i < n_hits; ++i) {
    if (math::abs(values[i] - mode) < tolerance ||
        values[i] == std::numeric_limits<scalar>::max()) {
      filtered_hits.push_back(hits[i]);
    }
  }
}

TRACCC_HOST_DEVICE
inline void sort_hits(const vecmem::device_vector<const scalar> values,
                      vecmem::device_vector<global_index_t> row) {
  assert(values.capacity() == row.capacity());
  const unsigned int n = row.capacity();

  for (unsigned int i = 1; i < n; ++i) {
    const unsigned int index = row[i];
    const scalar value = values[index];
    unsigned int j = i;

    while (j > 0 && value < values[row[j - 1]]) {
      row[j] = row[j - 1];
      --j;
    }
    row[j] = index;
  }
}

TRACCC_HOST_DEVICE
inline void make_seeds(
    const global_index_t globalIndex, const unsigned int n_images,
    const vecmem::data::vector_view<const scalar>& spacepoints_r_view,
    const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
    const vecmem::data::jagged_vector_view<const global_index_t>& hits_view,
    const vecmem::data::jagged_vector_view<const scalar>& cot_theta_matrix_view,
    const vecmem::data::vector_view<const scalar>& mode_values_view,
    const vecmem::data::vector_view<const unsigned int>& best_matches_view,
    vecmem::data::jagged_vector_view<global_index_t>& filtered_hits_view,
    edm::seed_collection::view& seeds_view) {
  if (globalIndex >= n_images) {
    return;
  }

  // does not pass the NN selection
  const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
  if (nn_outcome[globalIndex] == 0) {
    return;
  }

  const vecmem::jagged_device_vector<const global_index_t> hits(hits_view);
  const unsigned int n_hits = hits[globalIndex].capacity();
  if (n_hits < 3) {
    return;
  }

  // apply hit filtering
  const vecmem::device_vector<const scalar> spacepoints_r(spacepoints_r_view);
  const vecmem::jagged_device_vector<const scalar> cot_theta_matrix(
      cot_theta_matrix_view);
  const vecmem::device_vector<const scalar> mode_values(mode_values_view);
  const vecmem::device_vector<const unsigned int> best_matches(
      best_matches_view);
  vecmem::jagged_device_vector<global_index_t> filtered_hits(
      filtered_hits_view);

  const unsigned int n_filtered_hits = filtered_hits[globalIndex].capacity();
  if (n_filtered_hits < 3) {
    return;
  }

  const scalar mode = mode_values[globalIndex];
  const unsigned int best_match_index = best_matches[globalIndex];

  select_hits(&cot_theta_matrix[globalIndex][n_hits * best_match_index], mode,
              n_hits, hits[globalIndex], filtered_hits[globalIndex]);

  sort_hits(spacepoints_r, filtered_hits[globalIndex]);

  // create seeds
  edm::seed_collection::device seeds(seeds_view);

  seeds[globalIndex].bottom_index() = filtered_hits[globalIndex][0];
  seeds[globalIndex].middle_index() =
      filtered_hits[globalIndex][(n_filtered_hits - 1) / 2];
  seeds[globalIndex].top_index() =
      filtered_hits[globalIndex][n_filtered_hits - 1];
  seeds[globalIndex].quality() = 0.00f;
}

}  // namespace traccc::device
