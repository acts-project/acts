// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/definitions/math.hpp"
#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline unsigned long long getId(
    const vecmem::device_vector<const unsigned int>& shift_map,
    const unsigned int invalid_shift, const unsigned int index) {
  if (index >= shift_map.capacity())
    return 0ULL;
  const unsigned int shift = shift_map[index];
  if (shift == invalid_shift)
    return 0ULL;
  return 1ULL << shift;
}

template <typename bfield_t>
TRACCC_HOST_DEVICE inline void precompute_spacepoint_variables(
    const global_index_t globalIndex, const bfield_t& bfield,
    const unsigned int n_wedges,
    const vecmem::data::vector_view<const scalar>& wedges_view,
    const unsigned int invalid_shift,
    const vecmem::data::vector_view<const unsigned int>& shift_map_view,
    const edm::measurement_collection::const_view& measurements_view,
    const edm::spacepoint_collection::const_view::size_type n_spacepoints,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const scalar interactionCenter, const scalar interactionWindow,
    const scalar etaWidth, const unsigned long long module_veto,
    vecmem::data::vector_view<unsigned long long>&
        spacepoint_wedge_association_view,
    vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
    vecmem::data::vector_view<scalar>& spacepoint_r_view,
    vecmem::data::vector_view<scalar>& spacepoint_phi_view,
    vecmem::data::vector_view<scalar>& spacepoint_bfield_z_view) {
  if (globalIndex >= n_spacepoints) {
    return;
  }

  // defs
  const vecmem::device_vector<const scalar> wedges(wedges_view);
  const vecmem::device_vector<const unsigned int> shift_map(shift_map_view);
  // input
  const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
  // outputs
  vecmem::device_vector<unsigned long long> spacepoint_wedge_association(
      spacepoint_wedge_association_view);
  vecmem::device_vector<unsigned long long> spacepoint_detector_id(
      spacepoint_detector_id_view);
  vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
  vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);
  vecmem::device_vector<scalar> spacepoint_bfield_z(spacepoint_bfield_z_view);

  scalar z0_min = interactionCenter - interactionWindow;
  scalar z0_max = interactionCenter + interactionWindow;

  spacepoint_r[globalIndex] =
      hypotf(spacepoints[globalIndex].x(), spacepoints[globalIndex].y());
  scalar eta_at_min = math::asinh((spacepoints[globalIndex].z() - z0_min) /
                                  spacepoint_r[globalIndex]);
  scalar eta_at_max = math::asinh((spacepoints[globalIndex].z() - z0_max) /
                                  spacepoint_r[globalIndex]);

  spacepoint_phi[globalIndex] =
      atan2f(spacepoints[globalIndex].y(), spacepoints[globalIndex].x());
  spacepoint_phi[globalIndex] = ::remainderf(spacepoint_phi[globalIndex],
                                             2.0f * std::numbers::pi_v<scalar>);

  // Retrieve the corresponding measurement index and, from it, the geo id
  // all of this is to retrieve the volume_id and layer_id of the
  // space point
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const auto& meas =
      measurements.at(spacepoints[globalIndex].measurement_index_1());

  const auto volume_id = meas.surface_link().volume();
  const auto identifier_id = getId(shift_map, invalid_shift, volume_id);
  spacepoint_detector_id[globalIndex] = identifier_id;

  // veto layers
  if ((identifier_id & module_veto) != 0ULL) {
    return;
  }

  // understand what wedge(s) the space points belongs to
  // wedges define eta ranges and they may overlap
  //
  // this implementation has an upper limit on the number of wedges
  unsigned long long slice = 0ULL;
  for (unsigned int wedge_idx(0); wedge_idx < n_wedges; ++wedge_idx) {
    if (eta_at_min < wedges[wedge_idx] - etaWidth) {
      continue;
    }
    if (eta_at_max > wedges[wedge_idx] + etaWidth) {
      continue;
    }
    slice |= (1ULL << wedge_idx);
  }  // loop on wedges
  if (slice == 0ULL) {
    return;
  }
  spacepoint_wedge_association[globalIndex] = slice;

  // get the bfield
  const auto covfie_field_at_sp =
      bfield.at(spacepoints[globalIndex].x(), spacepoints[globalIndex].y(),
                spacepoints[globalIndex].z());
  spacepoint_bfield_z[globalIndex] = covfie_field_at_sp[2];
}
}  // namespace traccc::device
