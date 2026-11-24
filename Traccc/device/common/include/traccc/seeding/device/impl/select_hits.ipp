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
inline void sort_row(const vecmem::device_vector<const scalar> values,
                     vecmem::device_vector<unsigned int> row) {
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
inline unsigned int count_matches(const scalar* values,
                                  const unsigned int n_values,
                                  const scalar mode_value) {
  static constexpr scalar tolerance = 0.001f;
  unsigned int n_matches = 0u;
  for (unsigned int i = 0; i < n_values; ++i) {
    if (math::abs(values[i] - mode_value) > tolerance) {
      continue;
    }
    ++n_matches;
  }
  return n_matches;
}

TRACCC_HOST_DEVICE
inline void select_hits(
    const global_index_t globalIndex, const unsigned int n_images,
    const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
    const vecmem::data::jagged_vector_view<const unsigned int>& hits_view,
    const vecmem::data::jagged_vector_view<const scalar>& cot_theta_matrix_view,
    vecmem::data::jagged_vector_view<unsigned int>& indexes_view,
    vecmem::data::vector_view<scalar>& mode_values_view,
    vecmem::data::vector_view<unsigned int>& n_filtered_hits_view,
    vecmem::data::vector_view<unsigned int>& hit_indexes_view) {
  if (globalIndex >= n_images) {
    return;
  }

  // does not pass the NN selection
  const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
  if (nn_outcome[globalIndex] == 0u) {
    return;
  }

  const vecmem::jagged_device_vector<const unsigned int> hits(hits_view);
  const vecmem::jagged_device_vector<const scalar> cot_theta_matrix(
      cot_theta_matrix_view);
  vecmem::jagged_device_vector<unsigned int> indexes(indexes_view);
  vecmem::device_vector<unsigned int> hit_indexes(hit_indexes_view);
  vecmem::device_vector<unsigned int> n_filtered_hits(n_filtered_hits_view);
  vecmem::device_vector<scalar> mode_values(mode_values_view);

  // make the iota
  const unsigned int n_hits = indexes[globalIndex].capacity();
  for (unsigned int i = 0; i < n_hits; ++i) {
    indexes[globalIndex][i] = i;
  }

  // sort
  sort_row(cot_theta_matrix[globalIndex], indexes[globalIndex]);

  // get mode value
  scalar mode_value = std::numeric_limits<scalar>::max();
  unsigned int mode_counts = 0u;
  scalar previous_value = std::numeric_limits<scalar>::max();
  unsigned int current_counts = 0u;

  static constexpr scalar tolerance = 0.001f;

  const unsigned int n = indexes[globalIndex].capacity();
  for (unsigned int i = 0; i < n; ++i) {
    const unsigned int index = indexes[globalIndex][i];
    const scalar value = cot_theta_matrix[globalIndex][index];
    if (math::abs(value - previous_value) > tolerance) {
      if (current_counts > mode_counts) {
        mode_counts = current_counts;
        mode_value = previous_value;
      }
      previous_value = value;
      current_counts = 0u;
    }
    if (value == std::numeric_limits<scalar>::max()) {
      break;
    }
    ++current_counts;
  }
  mode_values[globalIndex] = mode_value;

  // get index of sp that matches the most
  unsigned int n_matches = 0u;
  unsigned int best_match = std::numeric_limits<unsigned int>::max();

  const unsigned int n_associated_hits = hits[globalIndex].capacity();
  for (unsigned int i = 0u; i < n_associated_hits; ++i) {
    unsigned int current_matches =
        count_matches(&cot_theta_matrix[globalIndex][i * n_associated_hits],
                      n_associated_hits, mode_value) +
        1;
    assert(current_matches <= n_associated_hits);
    if (current_matches <= n_matches) {
      continue;
    }
    n_matches = current_matches;
    best_match = i;
  }
  hit_indexes[globalIndex] = best_match;
  n_filtered_hits[globalIndex] = n_matches;
}

}  // namespace traccc::device
