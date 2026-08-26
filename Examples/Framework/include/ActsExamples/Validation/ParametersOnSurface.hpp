// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"

#include <optional>

namespace Acts {
class Surface;
}

namespace ActsExamples {

/// The track-state parameters to extract.
enum class TrackParameterType {
  /// Parameters before the measurement update
  Predicted,
  /// Parameters after the measurement update
  Filtered,
  /// Smoothed parameters
  Smoothed,
  /// Smoothed parameters with the state's own measurement removed
  Unbiased,
};

/// Compute the truth bound track parameters on a surface from the simulated
/// hits of a measurement.
///
/// Position, direction, and time are averaged over the simulated hits, the
/// momentum is taken from the first one, and charge and hypothesis from the
/// given particle.
///
/// @note The hits are averaged regardless of which particle produced them, so
///       a merged cluster gives a mixture. Same as `RootTrackParameterWriter`.
///
/// @param gctx The geometry context
/// @param surface The surface to express the truth parameters on
/// @param measurementIndex The index of the measurement on the surface
/// @param particle The truth particle
/// @param simHits The simulated hits container
/// @param measurementSimHitsMap Map from measurement index to simulated hits
/// @param logger A logger for messages
/// @return The parameters without covariance, or nullopt without truth hits
std::optional<Acts::BoundTrackParameters> truthParametersOnSurface(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    Index measurementIndex, const SimParticle& particle,
    const SimHitContainer& simHits,
    const MeasurementSimHitsMap& measurementSimHitsMap,
    const Acts::Logger& logger);

/// Extract the reconstructed bound track parameters on the reference surface
/// of a track state.
///
/// @param state The track state to extract the parameters from
/// @param parameterType Which parameters; if not set, the best available ones
///        (smoothed, filtered, or predicted). `Unbiased` is never picked
///        implicitly and has to be requested
/// @param hypothesis The particle hypothesis for the parameters
/// @return The parameters with covariance, or nullopt without a reference
///         surface or matching parameters
std::optional<Acts::BoundTrackParameters> recoParametersOnSurface(
    const ConstTrackStateProxy& state,
    std::optional<TrackParameterType> parameterType,
    const Acts::ParticleHypothesis& hypothesis);

}  // namespace ActsExamples
