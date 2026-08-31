// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
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

/// Whether parameters of the given type were obtained using the track state's
/// own measurement.
///
/// @param parameterType The parameter type
/// @return True for `Filtered` and `Smoothed`, false otherwise
bool parametersUseOwnMeasurement(TrackParameterType parameterType);

/// The residual of track state parameters against the state's own calibrated
/// measurement.
struct MeasurementResidual {
  /// The bound parameters the measurement constrains.
  Acts::VariableBoundSubspaceHelper subspace;
  /// `reco - measurement` in the full bound space, zero outside @c subspace.
  Acts::BoundVector residual;
  /// The covariance of @c residual, zero outside @c subspace.
  Acts::BoundMatrix covariance;
};

/// Compute the residual of track state parameters against the calibrated
/// measurement of that same state.
///
/// Both live in the local frame of the state's surface, so this needs neither
/// a geometry context nor truth information.
///
/// The residual covariance is @f$ V + HPH^T @f$ for parameters that do not
/// use the state's own measurement and @f$ V - HPH^T @f$ for those that do,
/// with @f$ V @f$ the measurement covariance, @f$ P @f$ the parameter
/// covariance, and @f$ H @f$ the projector. The subtracted form is not
/// positive definite in general, so the caller has to handle a non-positive
/// variance.
///
/// @param state The track state holding the measurement
/// @param parameters The reconstructed parameters on the state's surface
/// @param parameterType Which parameters @p parameters are, deciding the sign
///        of @f$ HPH^T @f$
/// @return The residual, or nullopt if the state carries no measurement
std::optional<MeasurementResidual> measurementResidual(
    const ConstTrackStateProxy& state,
    const Acts::BoundTrackParameters& parameters,
    TrackParameterType parameterType);

}  // namespace ActsExamples
