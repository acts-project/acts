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
/// traccc::device::scan_block_offsets function
struct scan_block_offsets_payload {
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
   * @brief View object to the block_offset vector
   */
  vecmem::data::vector_view<const int> block_offsets_view;

  /**
   * @brief View object to the scanned block_offset vector
   */
  vecmem::data::vector_view<int> scanned_block_offsets_view;
};

}  // namespace traccc::device
