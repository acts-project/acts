// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfUtils.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"

namespace Acts {
namespace detail {

auto extractMultiComponentState(const MultiTrajectory &traj,
                                const std::vector<size_t> &tips,
                                const std::map<size_t, ActsScalar> &weights,
                                StatesType type)
    -> MultiComponentBoundTrackParameters<SinglyCharged> {
  throw_assert(
      !tips.empty(),
      "need at least one component to extract trajectory of type " << type);

  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> cmps;
  std::shared_ptr<const Surface> surface;

  for (auto &tip : tips) {
    const auto proxy = traj.getTrackState(tip);

    throw_assert(weights.find(tip) != weights.end(),
                 "Could not find weight for idx " << tip);

    switch (type) {
      case StatesType::ePredicted:
        cmps.push_back(
            {weights.at(tip), proxy.predicted(), proxy.predictedCovariance()});
        break;
      case StatesType::eFiltered:
        cmps.push_back(
            {weights.at(tip), proxy.filtered(), proxy.filteredCovariance()});
        break;
      case StatesType::eSmoothed:
        cmps.push_back(
            {weights.at(tip), proxy.smoothed(), proxy.smoothedCovariance()});
    }

    if (!surface) {
      surface = proxy.referenceSurface().getSharedPtr();
    } else {
      throw_assert(
          surface->geometryId() == proxy.referenceSurface().geometryId(),
          "surface mismatch");
    }
  }

  return MultiComponentBoundTrackParameters<SinglyCharged>(surface, cmps);
}

void computePosteriorWeights(const MultiTrajectory &mt,
                             const std::vector<std::size_t> &tips,
                             std::map<std::size_t, double> &weights) {
  // Helper Function to compute detR
  auto computeDetR = [](const auto &trackState) -> ActsScalar {
    const auto predictedCovariance = trackState.predictedCovariance();

    return visit_measurement(
        trackState.calibrated(), trackState.calibratedCovariance(),
        trackState.calibratedSize(),
        [&](const auto calibrated, const auto calibratedCovariance) {
          constexpr size_t kMeasurementSize =
              decltype(calibrated)::RowsAtCompileTime;
          const auto H =
              trackState.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>()
                  .eval();

          return (H * predictedCovariance * H.transpose() +
                  calibratedCovariance)
              .determinant();
        });
  };

  // Find minChi2, this can be used to factor some things later in the
  // exponentiation
  const auto minChi2 =
      mt.getTrackState(*std::min_element(tips.begin(), tips.end(),
                                         [&](const auto &a, const auto &b) {
                                           return mt.getTrackState(a).chi2() <
                                                  mt.getTrackState(b).chi2();
                                         }))
          .chi2();

  // Loop over the tips and compute new weights
  for (auto tip : tips) {
    const auto state = mt.getTrackState(tip);
    const double chi2 = state.chi2() - minChi2;
    const double detR = computeDetR(state);

    // If something is not finite here, just leave the weight as it is
    if (std::isfinite(chi2) && std::isfinite(detR)) {
      const auto factor = std::sqrt(1. / detR) * std::exp(-0.5 * chi2);
      weights.at(tip) *= factor;
    }
  }
}

}  // namespace detail
}  // namespace Acts
