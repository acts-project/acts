// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include "Acts/TrackFitting/detail/SymmetricKlDistanceMatrix.hpp"

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceImpl(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
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
  std::sort(cmpCache.begin(), cmpCache.end(),
            [&](const auto &a, const auto &b) {
              return proj(a).weight < proj(b).weight;
            });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceAggressiveImpl(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const proj_t &proj,
    const angle_desc_t &desc) {
  Acts::detail::SymmetricKLDistanceMatrix distances(cmpCache, proj);

  auto remainingComponents = cmpCache.size();

  std::vector<std::tuple<double, std::size_t, std::size_t>> minDistPairs;
  minDistPairs.reserve(distances.size());
  std::vector<std::size_t> toSkip;
  while (remainingComponents > maxCmpsAfterMerge) {
    minDistPairs.clear();
    toSkip.clear();
    distances.minDistancePairs(minDistPairs);

    for (const auto &[_, minI, minJ] : minDistPairs) {
      // Check if we have touched one of these components already
      if (std::find(toSkip.begin(), toSkip.end(), minI) != toSkip.end() ||
          std::find(toSkip.begin(), toSkip.end(), minJ) != toSkip.end()) {
        continue;
      }

      // Set one component and compute associated distances
      cmpCache[minI] =
          mergeComponents(cmpCache[minI], cmpCache[minJ], proj, desc);
      distances.recomputeAssociatedDistances(minI, cmpCache, proj);

      // Since we modified this component in this pass, we shouldn't touch it
      // again
      toSkip.push_back(minI);

      // Set weight of the other component to -1 so we can remove it later and
      // mask its distances
      proj(cmpCache[minJ]).weight = -1.0;
      distances.maskAssociatedDistances(minJ);

      remainingComponents--;

      // Break if reached the required amount of components
      if (remainingComponents == maxCmpsAfterMerge) {
        break;
      }
    }
  }

  // Remove all components which are labeled with weight -1
  std::sort(cmpCache.begin(), cmpCache.end(),
            [&](const auto &a, const auto &b) {
              return proj(a).weight < proj(b).weight;
            });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

namespace Acts {

void reduceMixtureLargestWeights(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  std::sort(cmpCache.begin(), cmpCache.end(),
            [](const auto &a, const auto &b) { return a.weight > b.weight; });
  cmpCache.resize(maxCmpsAfterMerge);
}

void reduceMixtureWithKLDistance(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &surface) {
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

void reduceMixtureWithKLDistanceAggressive(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &surface) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  auto proj = [](auto &a) -> decltype(auto) { return a; };

  // We must differ between surface types, since there can be different
  // local coordinates
  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    reduceWithKLDistanceAggressiveImpl(cmpCache, maxCmpsAfterMerge, proj, desc);
  });
}

}  // namespace Acts
