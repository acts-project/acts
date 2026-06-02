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
#include "Acts/MagneticField/MagneticFieldContext.hpp"

#include <memory>
#include <optional>

namespace Acts {
class BasePropagator;
class PerigeeSurface;
}  // namespace Acts

namespace ActsExamples {

/// Express truth track parameters at a perigee surface.
///
/// Charged particles are propagated along their truth helix to @p surface
/// using @p propagator (via @c BasePropagator::propagateToSurface, which forces
/// the perigee surface to be reached). Neutral particles, which have no helix,
/// are linearly extrapolated to the straight-line intersection with @p surface;
/// for them the returned direction and momentum equal the input truth values.
///
/// @param geoContext Geometry context for the surface operations.
/// @param magFieldContext Magnetic field context for the propagation.
/// @param propagator Propagator used for charged particles.
/// @param surface Perigee surface to express the parameters on.
/// @param start Truth (curvilinear) parameters at the production vertex.
///
/// @return The bound parameters on @p surface, or @c std::nullopt if the
///         extrapolation (neutral) or propagation (charged) failed.
std::optional<Acts::BoundTrackParameters> propagateToPerigee(
    const Acts::GeometryContext& geoContext,
    const Acts::MagneticFieldContext& magFieldContext,
    const Acts::BasePropagator& propagator,
    std::shared_ptr<const Acts::PerigeeSurface> surface,
    const Acts::BoundTrackParameters& start);

}  // namespace ActsExamples
