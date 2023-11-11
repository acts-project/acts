// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include "Acts/TrackFitting/detail/MergeGaussianMixture.hpp"
#include "Acts/TrackFitting/detail/SymmetricKlDistanceMatrix.hpp"

template <typename component_t, typename component_projector_t,
          typename angle_desc_t>
auto mergeComponents(const component_t &a, const component_t &b,
                     const component_projector_t &proj,
                     const angle_desc_t &angle_desc) {
  assert(proj(a).weight >= 0.0 && proj(b).weight >= 0.0 &&
         "non-positive weight");

  std::array range = {std::ref(proj(a)), std::ref(proj(b))};
  const auto refProj = [](auto &c) {
    return std::tie(c.get().weight, c.get().boundPars, c.get().boundCov);
  };

  const auto mergedPars =
      Acts::detail::computeMixtureMean(range, refProj, angle_desc);
  const auto mergedCov = Acts::detail::computeMixtureCovariance(
      range, mergedPars, refProj, angle_desc);

  component_t ret = a;
  proj(ret).boundPars = mergedPars;
  proj(ret).boundCov = mergedCov;
  proj(ret).weight = proj(a).weight + proj(b).weight;

  return ret;
}

template <typename distance_t, typename proj_t, typename angle_desc_t>
void reduceWithDistanceImpl(std::vector<Acts::GsfComponent> &cmpCache,
                              std::size_t maxCmpsAfterMerge, const proj_t &proj,
                              const angle_desc_t &desc) {
  Acts::detail::SymmetricKLDistanceMatrix<distance_t>
      distances(cmpCache, proj);

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
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

namespace Acts {

void reduceMixtureLargestWeights(std::vector<GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  std::sort(cmpCache.begin(), cmpCache.end(),
            [](const auto &a, const auto &b) { return a.weight > b.weight; });
  cmpCache.resize(maxCmpsAfterMerge);
}

void reduceMixtureWithKLDistance(std::vector<GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &surface) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  auto proj = [](auto &a) -> decltype(auto) { return a; };

  // We must differ between surface types, since there can be different
  // local coordinates
  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    reduceWithDistanceImpl<Acts::detail::SymmetricKLDistanceQoP>(cmpCache, maxCmpsAfterMerge, proj, desc);
  });
}

}  // namespace Acts
