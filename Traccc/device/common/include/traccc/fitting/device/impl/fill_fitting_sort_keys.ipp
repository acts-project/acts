// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void fill_fitting_sort_keys(
    const global_index_t globalIndex,
    const edm::track_collection<default_algebra>::const_view&
        track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {
  const edm::track_collection<default_algebra>::const_device track_candidates(
      track_candidates_view);

  // Keys
  vecmem::device_vector<device::sort_key> keys_device(keys_view);

  // Param id
  vecmem::device_vector<unsigned int> ids_device(ids_view);

  if (globalIndex >= keys_device.size()) {
    return;
  }

  // Key = The number of measurements
  keys_device.at(globalIndex) = static_cast<traccc::scalar>(
      track_candidates.at(globalIndex).constituent_links().size());
  ids_device.at(globalIndex) = globalIndex;
}

}  // namespace traccc::device
