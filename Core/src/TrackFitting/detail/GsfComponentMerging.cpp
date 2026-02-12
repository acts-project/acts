// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"

namespace Acts {

std::tuple<BoundVector, BoundMatrix> detail::Gsf::mergeGaussianMixture(
    std::span<const GsfComponent> mixture, const Surface &surface,
    ComponentMergeMethod method) {
  return mergeGaussianMixture(mixture, std::identity{}, surface, method);
}

GsfComponent detail::Gsf::mergeTwoComponents(const GsfComponent &a,
                                             const GsfComponent &b,
                                             const Surface &surface) {
  assert(a.weight >= 0.0 && b.weight >= 0.0 && "non-positive weight");

  std::array components = {&a, &b};
  const auto proj = [](const GsfComponent *c) {
    return std::tie(c->weight, c->boundPars, c->boundCov);
  };
  auto [mergedPars, mergedCov] =
      angleDescriptionSwitch(surface, [&](const auto &desc) {
        return mergeGaussianMixtureMeanCov(components, proj, desc);
      });

  GsfComponent ret = a;
  ret.boundPars = mergedPars;
  ret.boundCov = mergedCov;
  ret.weight = a.weight + b.weight;
  return ret;
}

}  // namespace Acts
