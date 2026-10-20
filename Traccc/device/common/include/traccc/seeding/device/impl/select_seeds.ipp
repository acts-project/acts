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
#include "traccc/seeding/detail/triplet_sorter.hpp"
#include "traccc/seeding/seed_confirmation_helper.hpp"
#include "traccc/seeding/seed_selecting_helper.hpp"

// System include(s)
#include <cassert>
#include <cmath>
#include <limits>

namespace traccc::device {

namespace details {
// Compare two device triplets
struct device_triplet_comparator {
  constexpr bool operator()(const device_triplet& lhs,
                            const device_triplet& rhs) {
    if (lhs.weight != rhs.weight) {
      return lhs.weight > rhs.weight;
    } else {
      return std::tuple<unsigned int, unsigned int, unsigned int>(
                 lhs.spB, lhs.spM, lhs.spT) >
             std::tuple<unsigned int, unsigned int, unsigned int>(
                 rhs.spB, rhs.spM, rhs.spT);
    }
  }
};

// Finding minimum element algorithm
template <typename Comparator>
TRACCC_HOST_DEVICE std::size_t min_elem(const device_triplet* arr,
                                        const std::size_t begin_idx,
                                        const std::size_t end_idx,
                                        Comparator comp) {
  assert(begin_idx <= end_idx);
  std::size_t min_i = begin_idx;
  std::size_t next = begin_idx;
  for (std::size_t i = begin_idx + 1; i < end_idx; ++i) {
    ++next;
    if (comp(arr[min_i], arr[next])) {
      min_i = next;
    }
  }
  return min_i;
}

// Sorting algorithm for sorting seeds in the local memory
template <typename Comparator>
TRACCC_HOST_DEVICE void insertionSort(device_triplet* arr,
                                      const unsigned int begin_idx,
                                      const unsigned int n, Comparator comp) {
  int j = 0;
  device_triplet key = arr[begin_idx];
  for (unsigned int i = 0; i < n; ++i) {
    key = arr[begin_idx + i];
    j = static_cast<int>(i) - 1;
    while (j >= 0 &&
           !comp(arr[begin_idx + static_cast<unsigned int>(j)], key)) {
      arr[begin_idx + static_cast<unsigned int>(j + 1)] =
          arr[begin_idx + static_cast<unsigned int>(j)];
      j = j - 1;
    }
    arr[begin_idx + static_cast<unsigned int>(j + 1)] = key;
  }
}
/// Keep the @c max_n best triplets of a bounded array
TRACCC_HOST_DEVICE inline void push_bounded(device_triplet* arr,
                                            unsigned int& n,
                                            const unsigned int max_n,
                                            const device_triplet& t) {
  if (max_n == 0u) {
    return;
  }
  if (n < max_n) {
    arr[n++] = t;
    return;
  }
  const std::size_t min_i = min_elem(arr, 0, n, device_triplet_comparator{});
  if (device_triplet_comparator{}(t, arr[min_i])) {
    arr[min_i] = t;
  }
}

/// Seed selection of one middle spacepoint following the Acts triplet seed
/// filter with seed confirmation. The triplets of the middle spacepoint are
/// stored contiguously and grouped by their bottom spacepoint.
TRACCC_HOST_DEVICE inline void select_seeds_confirmed(
    const seedfinder_config& finder_config,
    const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_device& spacepoints,
    const triplet_counter_collection_types::const_device& triplet_counts,
    const device_triplet_collection_types::const_device& triplets,
    const edm::spacepoint_collection::const_device::const_proxy_type& spM,
    const unsigned int begin, const unsigned int end, device_triplet* data,
    edm::seed_collection::device& seeds_device) {
  // Candidates: quality seeds first, then the others
  device_triplet* quality = data;
  device_triplet* others = data + filter_config.maxQualitySeedsPerSpMConf;
  unsigned int n_quality = 0u;
  unsigned int n_others = 0u;

  unsigned int i = begin;
  while (i < end) {
    const device_triplet first = triplets[i];

    // Guard against the sentinel value for invalid triplets. The remaining
    // triplets of this group are skipped through its triplet counter below.
    if (first.spB == std::numeric_limits<unsigned int>::max()) [[unlikely]] {
      ++i;
      continue;
    }

    const unsigned int group_n =
        triplet_counts.at(static_cast<unsigned int>(first.counter_link))
            .m_nTriplets;
    const unsigned int group_end = i + group_n;
    const edm::spacepoint_collection::const_device::const_proxy_type spB =
        spacepoints.at(first.spB);
    const scalar spB_r = spB.radius();

    if (group_n < seed_confirmation_helper::min_top_candidates(
                      filter_config, spB_r, spM.z(), n_quality)) {
      i = group_end;
      continue;
    }
    const seed_confirmation_range& rng =
        seed_confirmation_helper::range(filter_config, spB.z());

    device_triplet best_other{};
    bool has_other = false;
    for (unsigned int j = i; j < group_end; ++j) {
      device_triplet t = triplets[j];
      if (!seed_confirmation_helper::confirm(
              filter_config, rng, spB_r, t.impact_parameter, t.z_vertex,
              t.n_compat, n_quality, t.weight, t.quality)) {
        continue;
      }
      if (t.quality) {
        push_bounded(quality, n_quality,
                     filter_config.maxQualitySeedsPerSpMConf, t);
      } else if (!has_other || t.weight > best_other.weight) {
        best_other = t;
        has_other = true;
      }
    }
    if (has_other && n_quality == 0u) {
      push_bounded(others, n_others, filter_config.maxSeedsPerSpMConf,
                   best_other);
    }
    i = group_end;
  }

  // Only quality seeds are kept if any were found
  device_triplet* selected = n_quality > 0u ? quality : others;
  const unsigned int n_selected = n_quality > 0u ? n_quality : n_others;
  insertionSort(selected, 0, n_selected, device_triplet_comparator{});
  const unsigned int max_seeds =
      math::min(n_selected, finder_config.maxSeedsPerSpM + 1u);
  for (unsigned int k = 0; k < max_seeds; ++k) {
    const device_triplet& t = selected[k];
    seeds_device.push_back({t.spB, t.spM, t.spT, static_cast<float>(t.weight)});
  }
}

}  // namespace details

// Select seeds kernel
TRACCC_HOST_DEVICE
inline void select_seeds(
    const global_index_t globalIndex, const seedfinder_config& finder_config,
    const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const traccc::details::spacepoint_grid_types::const_view& sp_view,
    const triplet_counter_spM_collection_types::const_view& spM_tc_view,
    const triplet_counter_collection_types::const_view& tc_view,
    const device_triplet_collection_types::const_view& triplet_view,
    device_triplet* data, edm::seed_collection::view seed_view) {
  // Check if anything needs to be done.
  const triplet_counter_spM_collection_types::const_device triplet_counts_spM(
      spM_tc_view);
  if (globalIndex >= triplet_counts_spM.size()) {
    return;
  }

  // Set up the device containers
  const triplet_counter_collection_types::const_device triplet_counts(tc_view);
  const edm::spacepoint_collection::const_device spacepoints{spacepoints_view};
  const traccc::details::spacepoint_grid_types::const_device sp_device(sp_view);

  const device_triplet_collection_types::const_device triplets(triplet_view);
  edm::seed_collection::device seeds_device(seed_view);

  // Current work item = middle spacepoint
  const triplet_counter_spM spM_counter = triplet_counts_spM.at(globalIndex);
  const sp_location spM_loc = spM_counter.spM;
  const unsigned int spM_idx = sp_device.bin(spM_loc.bin_idx)[spM_loc.sp_idx];
  const edm::spacepoint_collection::const_device::const_proxy_type spM =
      spacepoints.at(spM_idx);

  // Number of triplets added for this spM
  unsigned int n_triplets_per_spM = 0;

  const unsigned int end_triplets_spM =
      spM_counter.posTriplets + spM_counter.m_nTriplets;

  if (filter_config.seedConfirmation) {
    details::select_seeds_confirmed(
        finder_config, filter_config, spacepoints, triplet_counts, triplets,
        spM, spM_counter.posTriplets, end_triplets_spM, data, seeds_device);
    return;
  }

  // iterate over the triplets in the bin
  for (unsigned int i = spM_counter.posTriplets; i < end_triplets_spM; ++i) {
    device_triplet aTriplet = triplets[i];

    // spacepoints bottom for this triplet
    const unsigned int spB_idx = aTriplet.spB;

    // Guard against the sentinel value for invalid triplets.
    if (spB_idx == std::numeric_limits<unsigned int>::max()) [[unlikely]] {
      continue;
    }

    const edm::spacepoint_collection::const_device::const_proxy_type spB =
        spacepoints.at(spB_idx);
    // spacepoints top for this triplet
    const unsigned int spT_idx = aTriplet.spT;
    const edm::spacepoint_collection::const_device::const_proxy_type spT =
        spacepoints.at(spT_idx);

    // update weight of triplet and check if it is a good triplet
    if (filter_config.useTracccSeedCuts) {
      seed_selecting_helper::seed_weight(filter_config, spM, spB, spT,
                                         aTriplet.weight);
      if (!seed_selecting_helper::single_seed_cut(filter_config, spM, spB, spT,
                                                  aTriplet.weight)) {
        continue;
      }
    }

    // if the number of good triplets is larger than the threshold,
    // the triplet with the lowest weight is removed
    if (n_triplets_per_spM >= finder_config.maxSeedsPerSpM) {
      const std::size_t min_index =
          details::min_elem(data, 0, finder_config.maxSeedsPerSpM,
                            details::device_triplet_comparator{});

      if (details::device_triplet_comparator{}(aTriplet, data[min_index])) {
        data[min_index] = aTriplet;
      }
    }

    // if the number of good triplets is below the threshold, add
    // the current triplet to the array
    else if (n_triplets_per_spM < finder_config.maxSeedsPerSpM) {
      data[n_triplets_per_spM] = aTriplet;
      n_triplets_per_spM++;
    }
  }

  // sort the triplets per spM
  details::insertionSort(data, 0, n_triplets_per_spM,
                         details::device_triplet_comparator{});

  // the number of good seed per compatible middle spacepoint
  unsigned int n_seeds_per_spM = 0;

  // iterate over the good triplets for final selection of seeds
  for (unsigned int i = 0; i < n_triplets_per_spM; ++i) {
    const device_triplet& aTriplet = data[i];

    // if the number of seeds reaches the threshold, break
    if (n_seeds_per_spM >= finder_config.maxSeedsPerSpM + 1) {
      break;
    }

    // check if it is a good triplet
    if (!filter_config.useTracccSeedCuts ||
        seed_selecting_helper::cut_per_middle_sp(
            filter_config, spacepoints.at(aTriplet.spB), aTriplet.weight) ||
        n_seeds_per_spM == 0) {
      n_seeds_per_spM++;

      seeds_device.push_back({aTriplet.spB, aTriplet.spM, aTriplet.spT,
                              static_cast<float>(aTriplet.weight)});
    }
  }
}

}  // namespace traccc::device
