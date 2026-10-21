// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void fill_finding_propagation_sort_keys(
    const global_index_t globalIndex,
    const fill_finding_propagation_sort_keys_payload& payload) {
  const bound_track_parameters_collection_types::const_device params(
      payload.params_view);
  const vecmem::device_vector<const unsigned int> param_liveness(
      payload.param_liveness_view);

  // Keys
  vecmem::device_vector<device::sort_key> keys_device(payload.keys_view);

  // Param id
  vecmem::device_vector<unsigned int> ids_device(payload.ids_view);

  if (globalIndex >= keys_device.size()) {
    return;
  }

  /*
   * Adding a large constant factor to any dead tracks will ensure that they
   * all end up at the end of the array, and so they will produce minimal
   * thread divergence.
   */
  keys_device.at(globalIndex) =
      device::get_sort_key(params.at(globalIndex)) +
      (param_liveness.at(globalIndex) == 0u ? device::dead_track_sort_key_offset
                                            : 0.f);
  ids_device.at(globalIndex) = globalIndex;
}

}  // namespace traccc::device
