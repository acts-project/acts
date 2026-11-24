// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/seeding/seed_confirmation_helper.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void update_triplet_weights(
    const global_index_t globalIndex, const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const triplet_counter_spM_collection_types::const_view& spM_tc_view,
    const triplet_counter_collection_types::const_view& tc_view, scalar* data,
    device_triplet_collection_types::view triplet_view) {
  // Check if anything needs to be done.
  device_triplet_collection_types::device triplets(triplet_view);
  if (globalIndex >= triplets.size()) {
    return;
  }

  // Set up the device containers
  const edm::spacepoint_collection::const_device spacepoints{spacepoints_view};
  const triplet_counter_spM_collection_types::const_device triplet_counts_spM(
      spM_tc_view);
  const triplet_counter_collection_types::const_device triplet_counts(tc_view);

  // Current work item
  device_triplet this_triplet = triplets.at(globalIndex);
  const auto this_spT_idx = this_triplet.spT;

  if (this_spT_idx == std::numeric_limits<unsigned int>::max()) [[unlikely]] {
    return;
  }

  const edm::spacepoint_collection::const_device::const_proxy_type current_spM =
      spacepoints.at(this_triplet.spM);
  const edm::spacepoint_collection::const_device::const_proxy_type current_spT =
      spacepoints.at(this_spT_idx);

  const scalar currentTop_r = seed_confirmation_helper::top_radius(
      filter_config, current_spM.radius(), current_spM.z(),
      current_spT.radius(), current_spT.z());

  // if two compatible seeds with high distance in r are found, compatible
  // seeds span 5 layers
  // -> very good seed
  const scalar lowerLimitCurv =
      this_triplet.curvature - filter_config.deltaInvHelixDiameter;
  const scalar upperLimitCurv =
      this_triplet.curvature + filter_config.deltaInvHelixDiameter;
  std::size_t num_compat_seedR = 0;

  const triplet_counter mb_count =
      triplet_counts.at(static_cast<unsigned int>(this_triplet.counter_link));

  // Check if anything needs to be done.
  if (mb_count.m_nTriplets <= 1) {
    triplets.at(globalIndex).weight += seed_confirmation_helper::extra_weight(
        filter_config, 0u, this_triplet.impact_parameter, this_triplet.z_vertex,
        this_triplet.cot_theta);
    return;
  }

  // Get the position of the triplets which share this sabe midBottom doublet
  const unsigned int triplets_mb_begin =
      mb_count.posTriplets +
      triplet_counts_spM.at(mb_count.spM_counter_link).posTriplets;
  const unsigned int triplets_mb_end = triplets_mb_begin + mb_count.m_nTriplets;

  // iterate over triplets
  for (unsigned int i = triplets_mb_begin; i < triplets_mb_end; ++i) {
    // skip same triplet
    if (i == globalIndex) {
      continue;
    }

    const device_triplet other_triplet = triplets[i];

    const auto other_spT_idx = other_triplet.spT;

    if (other_spT_idx == std::numeric_limits<unsigned int>::max())
        [[unlikely]] {
      continue;
    }

    const edm::spacepoint_collection::const_device::const_proxy_type other_spT =
        spacepoints.at(other_spT_idx);

    // compared top SP should have at least deltaRMin distance
    const scalar otherTop_r = seed_confirmation_helper::top_radius(
        filter_config, current_spM.radius(), current_spM.z(),
        other_spT.radius(), other_spT.z());
    const scalar deltaR = currentTop_r - otherTop_r;
    if (math::fabs(deltaR) < filter_config.deltaRMin) {
      continue;
    }

    // curvature difference within limits?
    // TODO: how much slower than sorting all vectors by curvature
    // and breaking out of loop? i.e. is vector size large (e.g. in
    // jets?)
    if (other_triplet.curvature < lowerLimitCurv) {
      continue;
    }
    if (other_triplet.curvature > upperLimitCurv) {
      continue;
    }

    bool newCompSeed = true;

    for (std::size_t i_s = 0; i_s < num_compat_seedR; ++i_s) {
      const scalar previousDiameter = data[i_s];

      // original ATLAS code uses higher min distance for 2nd found
      // compatible seed (20mm instead of 5mm) add new compatible seed
      // only if distance larger than rmin to all other compatible
      // seeds
      if (math::fabs(previousDiameter - otherTop_r) < filter_config.deltaRMin) {
        newCompSeed = false;
        break;
      }
    }

    if (newCompSeed) {
      data[num_compat_seedR] = otherTop_r;
      this_triplet.weight += filter_config.compatSeedWeight;
      num_compat_seedR++;
    }

    if (num_compat_seedR >= filter_config.compatSeedLimit) {
      break;
    }
  }

  this_triplet.weight += seed_confirmation_helper::extra_weight(
      filter_config, static_cast<unsigned int>(num_compat_seedR),
      this_triplet.impact_parameter, this_triplet.z_vertex,
      this_triplet.cot_theta);
  triplets.at(globalIndex).weight = this_triplet.weight;
  triplets.at(globalIndex).n_compat =
      static_cast<unsigned int>(num_compat_seedR);
}

}  // namespace traccc::device
