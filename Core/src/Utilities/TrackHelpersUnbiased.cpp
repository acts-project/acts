// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

// This translation unit exists solely to instantiate the (very expensive)
// calculateUnbiasedParametersCovariance() once, against the type-erased
// AnyConstTrackStateProxy. Callers use the non-template overload declared in
// AnyTrackStateProxy.hpp and thus avoid re-instantiating the Eigen-heavy body
// for their own concrete track state proxy type.

namespace Acts {

std::pair<BoundVector, BoundMatrix> calculateUnbiasedParametersCovariance(
    const AnyConstTrackStateProxy& trackState) {
  // Delegates to the deprecated templated implementation. This is the one
  // place it is meant to be instantiated, so silence the deprecation warning
  // for the internal call only.
  // ACTS_PUSH_IGNORE_DEPRECATED()
  // return calculateUnbiasedParametersCovariance<AnyConstTrackStateProxy>(
  //     trackState);
  // ACTS_POP_IGNORE_DEPRECATED()

  if (!trackState.hasSmoothed()) {
    throw std::invalid_argument("track state has no smoothed parameters");
  }
  if (!trackState.hasCalibrated()) {
    throw std::invalid_argument("track state has no calibrated parameters");
  }

  return visit_measurement(
      trackState.calibratedSize(),
      [&]<std::size_t measdim>(std::integral_constant<std::size_t, measdim>) {
        FixedBoundSubspaceHelper<measdim> subspaceHelper =
            trackState.template projectorSubspaceHelper<measdim>();

        // TODO use subspace helper for projection instead
        auto H = subspaceHelper.projector();
        auto s = trackState.smoothed();
        auto C = trackState.smoothedCovariance();
        auto m = trackState.template calibrated<measdim>();
        auto V = trackState.template calibratedCovariance<measdim>();
        auto K =
            (C * H.transpose() * (H * C * H.transpose() - V).inverse()).eval();
        BoundVector unbiasedParamsVec = s + K * (m - H * s);
        BoundMatrix unbiasedParamsCov = C - K * H * C;
        return std::make_pair(unbiasedParamsVec, unbiasedParamsCov);
      });
}

}  // namespace Acts
