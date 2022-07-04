// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/GsfError.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/Utilities/detail/gaussian_mixture_helpers.hpp"

namespace Acts {

namespace detail {

template <typename Function>
class ScopeGuard {
  Function f;

 public:
  ScopeGuard(Function &&fun) : f(std::move(fun)) {}
  ~ScopeGuard() { f(); }
};

/// @brief Smoothing function, which takes two ranges of
/// MultiTrajectory-indices and the corresponding projectors.
template <typename component_iterator_t, typename fwd_projector_t = Identity,
          typename bwd_projector_t = Identity>
auto bayesianSmoothing(component_iterator_t fwdBegin,
                       component_iterator_t fwdEnd,
                       component_iterator_t bwdBegin,
                       component_iterator_t bwdEnd,
                       fwd_projector_t fwdProjector = fwd_projector_t{},
                       bwd_projector_t bwdProjector = bwd_projector_t{}) {
  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      smoothedState;

  using ResType = Result<decltype(smoothedState)>;

  for (auto fwd = fwdBegin; fwd != fwdEnd; ++fwd) {
    const auto &[weight_a, pars_a, cov_a] = fwdProjector(*fwd);
    throw_assert(cov_a, "for now we require a covariance here");

    for (auto bwd = bwdBegin; bwd != bwdEnd; ++bwd) {
      const auto &[weight_b, pars_b, cov_b] = bwdProjector(*bwd);
      throw_assert(cov_b, "for now we require a covariance here");

      const auto summedCov = *cov_a + *cov_b;

      const auto K = *cov_a * summedCov.inverse();
      const auto new_pars = (pars_a + K * (pars_b - pars_a)).eval();
      const auto new_cov = (K * *cov_b).eval();

      const auto diff = pars_a - pars_b;
      const ActsScalar exponent = diff.transpose() * summedCov.inverse() * diff;

      const auto new_weight = std::exp(-0.5 * exponent) * weight_a * weight_b;

      if (new_weight == 0) {
        return ResType(GsfError::SmoothingFailed);
      }

      smoothedState.push_back({new_weight, new_pars, new_cov});
    }
  }

  normalizeWeights(smoothedState, [](auto &tuple) -> decltype(auto) {
    return std::get<double>(tuple);
  });

  throw_assert(weightsAreNormalized(
                   smoothedState,
                   [](const auto &tuple) { return std::get<double>(tuple); }),
               "smoothed state not normalized");

  return ResType(smoothedState);
}

/// @brief Projector type which maps a MultiTrajectory-Index to a tuple of
/// [weight, parameters, covariance]. Therefore, it contains a MultiTrajectory
/// and for now a std::map for the weights
template <StatesType type, typename traj_t>
struct MultiTrajectoryProjector {
  const MultiTrajectory<traj_t> &mt;
  const std::map<MultiTrajectoryTraits::IndexType, double> &weights;

  auto operator()(MultiTrajectoryTraits::IndexType idx) const {
    const auto proxy = mt.getTrackState(idx);
    switch (type) {
      case StatesType::ePredicted:
        return std::make_tuple(weights.at(idx), proxy.predicted(),
                               std::optional{proxy.predictedCovariance()});
      case StatesType::eFiltered:
        return std::make_tuple(weights.at(idx), proxy.filtered(),
                               std::optional{proxy.filteredCovariance()});
      case StatesType::eSmoothed:
        return std::make_tuple(weights.at(idx), proxy.smoothed(),
                               std::optional{proxy.smoothedCovariance()});
    }
  }
};

/// @brief This function takes two MultiTrajectory objects and corresponding
/// index lists (one of the backward pass, one of the forward pass), combines
/// them, applies smoothing, and returns a new, single-component MultiTrajectory
/// @tparam ReturnSmootedStates If set to true, the function returns not only
/// combined MultiTrajectory, but also a std::vector contianing the
/// component-wise smoothed states
/// TODO this function does not handle outliers correctly at the moment I think
/// TODO change std::vector< size_t > to boost::small_vector for better
/// performance
template <typename traj_t, bool ReturnSmootedStates = false>
auto smoothAndCombineTrajectories(
    const MultiTrajectory<traj_t> &fwd,
    const std::vector<MultiTrajectoryTraits::IndexType> &fwdStartTips,
    const std::map<MultiTrajectoryTraits::IndexType, double> &fwdWeights,
    const MultiTrajectory<traj_t> &bwd,
    const std::vector<MultiTrajectoryTraits::IndexType> &bwdStartTips,
    const std::map<MultiTrajectoryTraits::IndexType, double> &bwdWeights,
    LoggerWrapper logger = getDummyLogger()) {
  // This vector gets only filled if ReturnSmootedStates is true
  std::vector<std::pair<const Surface *,
                        std::vector<std::tuple<double, BoundVector,
                                               std::optional<BoundSymMatrix>>>>>
      smoothedStates;

  // Use backward trajectory as basic trajectory, so that final trajectory is
  // ordered correctly. We ensure also that they are unique.
  std::vector<MultiTrajectoryTraits::IndexType> bwdTips = bwdStartTips;

  // Ensures that the bwd tips are unique and do not contain kInvalid which
  // represents an invalid trajectory state
  auto sortUniqueValidateBwdTips = [&]() {
    std::sort(bwdTips.begin(), bwdTips.end());
    bwdTips.erase(std::unique(bwdTips.begin(), bwdTips.end()), bwdTips.end());

    auto invalid_it = std::find(bwdTips.begin(), bwdTips.end(),
                                MultiTrajectoryTraits::kInvalid);
    if (invalid_it != bwdTips.end()) {
      bwdTips.erase(invalid_it);
    }
  };

  sortUniqueValidateBwdTips();

  KalmanFitterResult<traj_t> result;
  // result.fittedStates = std::make_shared<traj_t>();

  while (!bwdTips.empty()) {
    // Ensure that we update the bwd tips whenever we go to the next iteration
    // (This allows using continue etc in the loop)
    ScopeGuard scopeGuard([&]() {
      for (auto &tip : bwdTips) {
        const auto p = bwd.getTrackState(tip);
        tip = p.previous();
      }

      sortUniqueValidateBwdTips();
    });

    const auto firstBwdState = bwd.getTrackState(bwdTips.front());
    const auto &currentSurface = firstBwdState.referenceSurface();

    // Search corresponding forward tips
    const auto bwdGeoId = currentSurface.geometryId();
    std::vector<MultiTrajectoryTraits::IndexType> fwdTips;

    for (const auto tip : fwdStartTips) {
      fwd.visitBackwards(tip, [&](const auto &state) {
        if (state.referenceSurface().geometryId() == bwdGeoId) {
          fwdTips.push_back(state.index());
        }
      });
    }

    // Check if we have forward tips
    if (fwdTips.empty()) {
      ACTS_WARNING("Did not find forward states on surface " << bwdGeoId);
      continue;
    }

    // Ensure we have no duplicates
    std::sort(fwdTips.begin(), fwdTips.end());
    fwdTips.erase(std::unique(fwdTips.begin(), fwdTips.end()), fwdTips.end());

    // Add state to MultiTrajectory
    result.lastTrackIndex = result.fittedStates.addTrackState(
        TrackStatePropMask::All, result.lastTrackIndex);
    result.processedStates++;

    auto proxy = result.fittedStates.getTrackState(result.lastTrackIndex);

    // This way we copy all relevant flags and the calibrated field. However
    // this assumes that the relevant flags do not differ between components
    proxy.copyFrom(firstBwdState);

    // Define some Projector types we need in the following
    using PredProjector =
        MultiTrajectoryProjector<StatesType::ePredicted, traj_t>;
    using FiltProjector =
        MultiTrajectoryProjector<StatesType::eFiltered, traj_t>;

    if (proxy.typeFlags().test(Acts::TrackStateFlag::HoleFlag)) {
      result.measurementHoles++;
    } else {
      // We also need to save outlier states here, otherwise they would not be
      // included in the MT if they are at the end of the track
      result.lastMeasurementIndex = result.lastTrackIndex;
    }

    // If we have a hole or an outlier, just take the combination of filtered
    // and predicted and no smoothed state
    if (not proxy.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
      const auto [mean, cov] = combineBoundGaussianMixture(
          bwdTips.begin(), bwdTips.end(), FiltProjector{bwd, bwdWeights});

      proxy.predicted() = mean;
      proxy.predictedCovariance() = cov.value();
      using PM = TrackStatePropMask;
      proxy.shareFrom(proxy, PM::Predicted, PM::Filtered);

    }
    // If we have a measurement, do the smoothing
    else {
      result.measurementStates++;

      // The predicted state is the forward pass
      const auto [fwdMeanPred, fwdCovPred] = combineBoundGaussianMixture(
          fwdTips.begin(), fwdTips.end(), PredProjector{fwd, fwdWeights});
      proxy.predicted() = fwdMeanPred;
      proxy.predictedCovariance() = fwdCovPred.value();

      // The filtered state is the backward pass
      const auto [bwdMeanFilt, bwdCovFilt] = combineBoundGaussianMixture(
          bwdTips.begin(), bwdTips.end(), FiltProjector{bwd, bwdWeights});
      proxy.filtered() = bwdMeanFilt;
      proxy.filteredCovariance() = bwdCovFilt.value();

      // Do the smoothing
      auto smoothedStateResult = bayesianSmoothing(
          fwdTips.begin(), fwdTips.end(), bwdTips.begin(), bwdTips.end(),
          PredProjector{fwd, fwdWeights}, FiltProjector{bwd, bwdWeights});

      if (!smoothedStateResult.ok()) {
        ACTS_WARNING("Smoothing failed on " << bwdGeoId);
        continue;
      }

      const auto &smoothedState = *smoothedStateResult;

      if constexpr (ReturnSmootedStates) {
        smoothedStates.push_back({&currentSurface, smoothedState});
      }

      // The smoothed state is a combination
      const auto [smoothedMean, smoothedCov] = combineBoundGaussianMixture(
          smoothedState.begin(), smoothedState.end());
      proxy.smoothed() = smoothedMean;
      proxy.smoothedCovariance() = smoothedCov.value();
      ACTS_VERBOSE("Added smoothed state to MultiTrajectory");
    }
  }

  result.smoothed = true;
  result.reversed = true;
  result.finished = true;

  if constexpr (ReturnSmootedStates) {
    return std::make_tuple(result, smoothedStates);
  } else {
    return std::make_tuple(result);
  }
}

}  // namespace detail
}  // namespace Acts
