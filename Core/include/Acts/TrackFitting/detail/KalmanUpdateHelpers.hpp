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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts::detail {

/// This function encapsulates the Kalman update performed on a MultiTrajectory
/// for a single source link.
/// @tparam propagator_state_t The propagator state type
/// @tparam stepper_t The stepper type
/// @tparam extensions_t The type of the extensions used for the update
/// @param state The propagator state
/// @param stepper The stepper
/// @param extensions The extension used for the update
/// @param surface The current surface
/// @param source_link The source-link used for the update
/// @param fittedStates The Multitrajectory to that we add the state
/// @param lastTrackIndex The parent index for the new state in the MT
/// @param doCovTransport Whether to perform a covariance transport when
/// computing the bound state or not
/// @param freeToBoundCorrection Correction for non-linearity effect during transform from free to bound (only corrected when performing CovTransport)
template <typename propagator_state_t, typename stepper_t,
          typename extensions_t, typename traj_t>
auto kalmanHandleMeasurement(
    const CalibrationContext &calibrationContext, propagator_state_t &state,
    const stepper_t &stepper, const extensions_t &extensions,
    const Surface &surface, const SourceLink &source_link, traj_t &fittedStates,
    const std::size_t lastTrackIndex, bool doCovTransport, const Logger &logger,
    const FreeToBoundCorrection &freeToBoundCorrection = FreeToBoundCorrection(
        false)) -> Result<typename traj_t::TrackStateProxy> {
  // Add a <mask> TrackState entry multi trajectory. This allocates storage for
  // all components, which we will set later.
  TrackStatePropMask mask =
      TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
      TrackStatePropMask::Smoothed | TrackStatePropMask::Jacobian |
      TrackStatePropMask::Calibrated;
  typename traj_t::TrackStateProxy trackStateProxy =
      fittedStates.makeTrackState(mask, lastTrackIndex);

  // Set the trackStateProxy components with the state from the ongoing
  // propagation
  {
    trackStateProxy.setReferenceSurface(surface.getSharedPtr());
    // Bind the transported state to the current surface
    auto res = stepper.boundState(state.stepping, surface, doCovTransport,
                                  freeToBoundCorrection);
    if (!res.ok()) {
      ACTS_ERROR("Propagate to surface " << surface.geometryId()
                                         << " failed: " << res.error());
      return res.error();
    }
    const auto &[boundParams, jacobian, pathLength] = *res;

    // Fill the track state
    trackStateProxy.predicted() = boundParams.parameters();
    trackStateProxy.predictedCovariance() = state.stepping.cov;

    trackStateProxy.jacobian() = jacobian;
    trackStateProxy.pathLength() = pathLength;
  }

  // We have predicted parameters, so calibrate the uncalibrated input
  // measurement
  extensions.calibrator(state.geoContext, calibrationContext, source_link,
                        trackStateProxy);

  // Get and set the type flags
  {
    auto typeFlags = trackStateProxy.typeFlags();
    typeFlags.set(TrackStateFlag::ParameterFlag);
    if (surface.surfaceMaterial() != nullptr) {
      typeFlags.set(TrackStateFlag::MaterialFlag);
    }

    // Check if the state is an outlier.
    // If not:
    // - run Kalman update
    // - tag it as a measurement
    // - update the stepping state.
    // Else, just tag it as an outlier
    if (!extensions.outlierFinder(trackStateProxy)) {
      // Run Kalman update
      auto updateRes = extensions.updater(state.geoContext, trackStateProxy,
                                          state.options.direction, logger);
      if (!updateRes.ok()) {
        ACTS_ERROR("Update step failed: " << updateRes.error());
        return updateRes.error();
      }
      // Set the measurement type flag
      typeFlags.set(TrackStateFlag::MeasurementFlag);
    } else {
      ACTS_VERBOSE(
          "Filtering step successful. But measurement is determined "
          "to be an outlier. Stepping state is not updated.")
      // Set the outlier type flag
      typeFlags.set(TrackStateFlag::OutlierFlag);
      trackStateProxy.shareFrom(trackStateProxy, TrackStatePropMask::Predicted,
                                TrackStatePropMask::Filtered);
    }
  }

  return trackStateProxy;
}

/// This function encapsulates what actions should be performed on a
/// MultiTrajectory when we have no measurement.
/// If there are no source links on surface, add either a hole or passive
/// material TrackState entry multi trajectory. No storage allocation for
/// uncalibrated/calibrated measurement and filtered parameter
/// @tparam propagator_state_t The propagator state type
/// @tparam stepper_t The stepper type
/// @param state The propagator state
/// @param stepper The stepper
/// @param surface The current surface
/// @param fittedStates The Multitrajectory to that we add the state
/// @param lastTrackIndex The parent index for the new state in the MT
/// @param doCovTransport Whether to perform a covariance transport when
/// computing the bound state or not
/// @param freeToBoundCorrection Correction for non-linearity effect during transform from free to bound (only corrected when performing CovTransport)
template <typename propagator_state_t, typename stepper_t, typename traj_t>
auto kalmanHandleNoMeasurement(
    propagator_state_t &state, const stepper_t &stepper, const Surface &surface,
    traj_t &fittedStates, const std::size_t lastTrackIndex, bool doCovTransport,
    const Logger &logger,
    const FreeToBoundCorrection &freeToBoundCorrection = FreeToBoundCorrection(
        false)) -> Result<typename traj_t::TrackStateProxy> {
  // Add a <mask> TrackState entry multi trajectory. This allocates storage for
  // all components, which we will set later.
  TrackStatePropMask mask = TrackStatePropMask::Predicted |
                            TrackStatePropMask::Smoothed |
                            TrackStatePropMask::Jacobian;
  typename traj_t::TrackStateProxy trackStateProxy =
      fittedStates.makeTrackState(mask, lastTrackIndex);

  // Set the trackStateProxy components with the state from the ongoing
  // propagation
  {
    trackStateProxy.setReferenceSurface(surface.getSharedPtr());
    // Bind the transported state to the current surface
    auto res = stepper.boundState(state.stepping, surface, doCovTransport,
                                  freeToBoundCorrection);
    if (!res.ok()) {
      return res.error();
    }
    const auto &[boundParams, jacobian, pathLength] = *res;

    // Fill the track state
    trackStateProxy.predicted() = boundParams.parameters();
    trackStateProxy.predictedCovariance() = state.stepping.cov;

    trackStateProxy.jacobian() = jacobian;
    trackStateProxy.pathLength() = pathLength;

    // Set the filtered parameter index to be the same with predicted
    // parameter
    trackStateProxy.shareFrom(trackStateProxy, TrackStatePropMask::Predicted,
                              TrackStatePropMask::Filtered);
  }

  // Set the track state flags
  auto typeFlags = trackStateProxy.typeFlags();
  typeFlags.set(TrackStateFlag::ParameterFlag);
  if (surface.surfaceMaterial() != nullptr) {
    typeFlags.set(TrackStateFlag::MaterialFlag);
  }
  if (surface.associatedDetectorElement() != nullptr) {
    ACTS_VERBOSE("Detected hole on " << surface.geometryId());
    // If the surface is sensitive, set the hole type flag
    typeFlags.set(TrackStateFlag::HoleFlag);
  } else if (surface.surfaceMaterial() != nullptr) {
    ACTS_VERBOSE("Detected in-sensitive surface " << surface.geometryId());
  }

  return trackStateProxy;
}

}  // namespace Acts::detail
