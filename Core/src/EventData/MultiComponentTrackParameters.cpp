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

#include <cstddef>
#include <ranges>

namespace Acts {

BoundTrackParameters MultiComponentBoundTrackParameters::merge(
    const ComponentMergeMethod method) const {
  assert(!hasCovariance() &&
         "Merging is only supported if covariance matrices are available");
  const auto [singleParams, singleCov] = detail::Gsf::mergeGaussianMixture(
      std::views::iota(std::size_t{0}, size()),
      [this](const std::size_t i)
          -> std::tuple<double, const BoundVector&, const BoundMatrix&> {
        return {m_weights[i], m_parameters[i], m_covariances[i]};
      },
      *m_surface, method);
  return BoundTrackParameters(m_surface, singleParams, singleCov,
                              m_particleHypothesis);
}

}  // namespace Acts
