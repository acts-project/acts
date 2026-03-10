// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/MultiComponentTrackParameters.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"

namespace Acts {

BoundTrackParameters MultiComponentBoundTrackParameters::merge(
    const ComponentMergeMethod method) const {
  const bool hasCov = std::get<2>(m_components.front()).has_value();

  if (!hasCov) {
    const BoundVector singleParams = detail::Gsf::mergeGaussianMixtureParams(
        m_components,
        [](const auto& cmp) -> std::tuple<double, const BoundVector&> {
          return {std::get<0>(cmp), std::get<1>(cmp)};
        },
        *m_surface, method);
    return BoundTrackParameters(m_surface, singleParams, std::nullopt,
                                m_particleHypothesis);
  }

  const auto [singleParams, singleCov] = detail::Gsf::mergeGaussianMixture(
      m_components,
      [](const auto& cmp)
          -> std::tuple<double, const BoundVector&, const BoundMatrix&> {
        return {std::get<0>(cmp), std::get<1>(cmp), std::get<2>(cmp).value()};
      },
      *m_surface, method);
  return BoundTrackParameters(m_surface, singleParams, singleCov,
                              m_particleHypothesis);
}

}  // namespace Acts
