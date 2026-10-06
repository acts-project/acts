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
/// traccc::device::fill_inverted_ids function
struct fill_inverted_ids_payload {
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
   * @brief View object to the inverted ids
   */
  vecmem::data::vector_view<unsigned int> inverted_ids_view;
};

}  // namespace traccc::device
