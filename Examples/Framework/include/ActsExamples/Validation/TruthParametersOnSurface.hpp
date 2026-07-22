// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"

#include <optional>

namespace Acts {
class Surface;
}

namespace ActsExamples {

/// Truth bound track parameters expressed on a specific surface.
struct TruthParametersOnSurface {
  /// The truth bound parameters loc0, loc1, phi, theta, q/p, t
  Acts::BoundVector params = Acts::BoundVector::Zero();
  /// The truth pseudorapidity
  double eta = 0;
  /// The truth azimuthal angle
  double phi = 0;
  /// The truth transverse momentum
  double pt = 0;
  /// The truth charge
  int charge = 0;
};

/// Compute the truth bound track parameters on a surface from the simulated
/// hits associated to a measurement.
///
/// The local position, direction, and time are obtained by averaging the
/// simulated hits on the surface. The momentum is taken from the first
/// simulated hit; the charge is provided by the caller, e.g. from the truth
/// particle matched to the track.
///
/// @param gctx The geometry context
/// @param surface The surface to express the truth parameters on
/// @param measurementIndex The index of the measurement on the surface
/// @param truthCharge The charge of the truth particle
/// @param simHits The simulated hits container
/// @param measurementSimHitsMap Map from measurement index to simulated hits
/// @param logger A logger for messages
/// @return The truth parameters, or nullopt if no simulated hits are
///         associated to the measurement
std::optional<TruthParametersOnSurface> truthParametersOnSurface(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    Index measurementIndex, double truthCharge, const SimHitContainer& simHits,
    const MeasurementSimHitsMap& measurementSimHitsMap,
    const Acts::Logger& logger);

}  // namespace ActsExamples
