// This file is part of the Acts project.
//
// Copyright (C) 2021-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
namespace detail {
/// This function creates and returns a track state proxy
/// @tparam propagator_state_t The propagator state type
/// @tparam stepper_t The stepper type
/// @param state The propagator state
/// @param stepper The stepper
/// @param surface The current surface
/// @param fittedStates The Multitrajectory to that we add the state
/// @param lastTrackIndex The parent index for the new state in the MT
/// @param doCovTransport Whether to perform a covariance transport when computing the bound state or not
/// @param freeToBoundCorrection Correction for non-linearity effect during transform from free to bound (only corrected when performing CovTransport)
/// @param getOnLastTrackIndex If set, no new TrackState will be added
template <typename propagator_state_t, typename stepper_t, typename traj_t>
auto getTrackStateProxy(propagator_state_t &state, const stepper_t &stepper,
                        const Surface &surface, traj_t &fittedStates,
                        const size_t lastTrackIndex, bool doCovTransport,
                        const Logger &logger,
                        const FreeToBoundCorrection &freeToBoundCorrection =
                            FreeToBoundCorrection(false),
                        TrackStatePropMask mask = TrackStatePropMask::All,
                        bool getOnLastTrackIndex = false)
    -> Result<typename traj_t::TrackStateProxy> {
  size_t currentTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;

  if (getOnLastTrackIndex) {
    if (lastTrackIndex < fittedStates.size() - 1 ||
        lastTrackIndex == Acts::MultiTrajectoryTraits::kInvalid) {
      ACTS_WARNING("No TrackState with index "
                   << lastTrackIndex << " found. Add new TrackState.");
      currentTrackIndex = fittedStates.addTrackState(mask, lastTrackIndex);
    }
    currentTrackIndex = lastTrackIndex;
  } else {
    // add a full TrackState entry multi trajectory
    // (this allocates storage for all components, we will set them later)
    currentTrackIndex = fittedStates.addTrackState(mask, lastTrackIndex);
  }

  // now get track state proxy back
  typename traj_t::TrackStateProxy trackStateProxy =
      fittedStates.getTrackState(currentTrackIndex);

  trackStateProxy.setReferenceSurface(surface.getSharedPtr());

  // Bind the transported state to the current surface
  auto res = stepper.boundState(state.stepping, surface, doCovTransport,
                                freeToBoundCorrection);
  if (!res.ok()) {
    ACTS_ERROR("Propagate to surface " << surface.geometryId()
                                       << " failed: " << res.error());
    return res.error();
  }
  auto &[boundParams, jacobian, pathLength] = *res;

  // Fill the track state
  trackStateProxy.predicted() = std::move(boundParams.parameters());
  if (boundParams.covariance().has_value()) {
    trackStateProxy.predictedCovariance() =
        std::move(*boundParams.covariance());
  }
  trackStateProxy.jacobian() = std::move(jacobian);
  trackStateProxy.pathLength() = std::move(pathLength);

  return trackStateProxy;
}

}  // namespace detail
}  // namespace Acts
