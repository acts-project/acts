// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/MultiComponentTrackParameters.hpp"

#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"

#include <cstddef>
#include <numeric>
#include <ranges>

namespace Acts {

BoundTrackParameters MultiComponentBoundTrackParameters::merge(
    const ComponentMergeMethod method) const {
  if (empty()) {
    throw std::logic_error(
        "Cannot merge MultiComponentBoundTrackParameters with zero components");
  }

  if (size() == 1) {
    return BoundTrackParameters(
        m_surface, m_parameters[0],
        hasCovariance() ? std::optional(m_covariances[0]) : std::nullopt,
        m_particleHypothesis);
  }

  if (!hasCovariance()) {
    const BoundVector singleParams = detail::Gsf::mergeGaussianMixtureParams(
        std::views::iota(std::size_t{0}, size()),
        [this](const std::size_t i) -> std::tuple<double, const BoundVector&> {
          return {m_weights[i], m_parameters[i]};
        },
        *m_surface, method);
    return BoundTrackParameters(m_surface, singleParams, std::nullopt,
                                m_particleHypothesis);
  }

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

void MultiComponentBoundTrackParameters::reserve(const std::size_t n) {
  m_weights.reserve(n);
  m_parameters.reserve(n);
  if (hasCovariance()) {
    m_covariances.reserve(n);
  }
}

void MultiComponentBoundTrackParameters::clear() {
  m_weights.clear();
  m_parameters.clear();
  m_covariances.clear();
}

void MultiComponentBoundTrackParameters::pushComponent(
    const double weight, const ParametersVector& params) {
  if (hasCovariance()) {
    throw std::logic_error(
        "Cannot push component without covariance to "
        "MultiComponentBoundTrackParameters with covariance");
  }

  m_weights.push_back(weight);
  m_parameters.push_back(params);
}

void MultiComponentBoundTrackParameters::pushComponent(
    const double weight, const ParametersVector& params,
    const CovarianceMatrix& cov) {
  if (!hasCovariance()) {
    throw std::logic_error(
        "Cannot push component with covariance to "
        "MultiComponentBoundTrackParameters without covariance");
  }

  m_weights.push_back(weight);
  m_parameters.push_back(params);
  m_covariances.push_back(cov);
}

void MultiComponentBoundTrackParameters::pushComponent(
    const double weight, const ParametersVector& params,
    const std::optional<CovarianceMatrix>& cov) {
  if (hasCovariance() != cov.has_value()) {
  }

  m_weights.push_back(weight);
  m_parameters.push_back(params);
  if (hasCovariance()) {
    m_covariances.push_back(*cov);
  }
}

void MultiComponentBoundTrackParameters::normalizeWeights() {
  const double sumWeights =
      std::accumulate(m_weights.begin(), m_weights.end(), 0.0);
  if (sumWeights <= 0.0) {
    throw std::logic_error(
        "Cannot normalize weights of MultiComponentBoundTrackParameters: sum "
        "of weights is not strictly positive: " +
        std::to_string(sumWeights));
  }
  for (double& weight : m_weights) {
    weight /= sumWeights;
  }
}

}  // namespace Acts
