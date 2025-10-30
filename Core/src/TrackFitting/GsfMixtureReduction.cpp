// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include "Acts/TrackFitting/detail/SymmetricKlDistanceMatrix.hpp"

#include <algorithm>

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceImpl(std::vector<Acts::GsfComponent> &cmpCache,
                              std::size_t maxCmpsAfterMerge, const proj_t &proj,
                              const angle_desc_t &desc) {
  Acts::detail::SymmetricKLDistanceMatrix distances(cmpCache, proj);

  auto remainingComponents = cmpCache.size();

  while (remainingComponents > maxCmpsAfterMerge) {
    const auto [minI, minJ] = distances.minDistancePair();

    // Set one component and compute associated distances
    cmpCache[minI] =
        mergeComponents(cmpCache[minI], cmpCache[minJ], proj, desc);
    distances.recomputeAssociatedDistances(minI, cmpCache, proj);

    // Set weight of the other component to -1 so we can remove it later and
    // mask its distances
    proj(cmpCache[minJ]).weight = -1.0;
    distances.maskAssociatedDistances(minJ);

    remainingComponents--;
  }

  // Remove all components which are labeled with weight -1
  std::ranges::sort(cmpCache, {},
                    [&](const auto &c) { return proj(c).weight; });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

namespace Acts {

void reduceMixtureLargestWeights(std::vector<GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface & /*surface*/) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  std::nth_element(
      cmpCache.begin(), cmpCache.begin() + maxCmpsAfterMerge, cmpCache.end(),
      [](const auto &a, const auto &b) { return a.weight > b.weight; });
  cmpCache.resize(maxCmpsAfterMerge);
}

void reduceMixtureWithKLDistance(std::vector<Acts::GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &surface) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  auto proj = [](auto &a) -> decltype(auto) { return a; };

  // We must differ between surface types, since there can be different
  // local coordinates
  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    reduceWithKLDistanceImpl(cmpCache, maxCmpsAfterMerge, proj, desc);
  });
}

void reduceMixtureWithKLDistanceNaive(std::vector<Acts::GsfComponent> &cmpCache,
                                      std::size_t maxCmpsAfterMerge,
                                      const Surface &surface) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  auto proj = [](auto &a) -> decltype(auto) { return a; };

  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    while (cmpCache.size() > maxCmpsAfterMerge) {
      // Recompute ALL distances every iteration (naive approach)
      double minDistance = std::numeric_limits<double>::max();
      std::size_t minI = 0;
      std::size_t minJ = 0;

      for (std::size_t i = 0; i < cmpCache.size(); ++i) {
        for (std::size_t j = i + 1; j < cmpCache.size(); ++j) {
          double distance = detail::computeSymmetricKlDivergence(
              cmpCache[i], cmpCache[j], proj);
          if (distance < minDistance) {
            minDistance = distance;
            minI = i;
            minJ = j;
          }
        }
      }

      // Merge the two closest components
      cmpCache[minI] = detail::mergeComponents(cmpCache[minI], cmpCache[minJ], proj, desc);

      // Remove the merged component immediately
      cmpCache.erase(cmpCache.begin() + minJ);
    }
  });

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

}  // namespace Acts
