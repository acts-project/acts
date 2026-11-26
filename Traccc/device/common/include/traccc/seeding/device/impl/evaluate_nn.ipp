// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void evaluate_nn(
    const global_index_t globalIndex, const unsigned int n_images,
    const vecmem::data::jagged_vector_view<const scalar>& images_view,
    const vecmem::data::vector_view<const scalar>& maxima_x_view,
    vecmem::data::vector_view<unsigned int>& nn_outcome_view) {
  if (globalIndex >= n_images) {
    return;
  }

  const vecmem::jagged_device_vector<const scalar> images(images_view);
  const vecmem::device_vector<const scalar> maxima_x(maxima_x_view);
  vecmem::device_vector<unsigned int> nn_outcome(nn_outcome_view);

  // For now we pass everything
  nn_outcome[globalIndex] = 1u;
}
}  // namespace traccc::device
