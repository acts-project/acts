// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
namespace traccc::device {

/// (Event Data) Payload for the @c
/// traccc::device::block_inclusive_scan function
struct block_inclusive_scan_payload {
  /**
   * @brief View object to the sorted track
   */
  vecmem::data::vector_view<const unsigned int> sorted_ids_view;

  /**
   * @brief Whether to terminate the calculation
   */
  int* terminate;

  /**
   * @brief The number of accepted tracks
   */
  unsigned int* n_accepted;

  /**
   * @brief The number of updated tracks
   */
  unsigned int* n_updated_tracks;

  /**
   * @brief View object to the whether track id is updated
   */
  vecmem::data::vector_view<const int> is_updated_view;

  /**
   * @brief View object to the block offset vector
   */
  vecmem::data::vector_view<int> block_offsets_view;

  /**
   * @brief View object to the prefix_sum vector
   */
  vecmem::data::vector_view<int> prefix_sums_view;
};

}  // namespace traccc::device
