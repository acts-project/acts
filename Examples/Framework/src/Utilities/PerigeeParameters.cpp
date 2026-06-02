// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/PerigeeParameters.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <utility>

namespace ActsExamples {

std::optional<Acts::BoundTrackParameters> propagateToPerigee(
    const Acts::GeometryContext& geoContext,
    const Acts::MagneticFieldContext& magFieldContext,
    const Acts::BasePropagator& propagator,
    std::shared_ptr<const Acts::PerigeeSurface> surface,
    const Acts::BoundTrackParameters& start) {
  const Acts::Vector3 direction = start.direction();

  // Closest intersection of the (straight) line through the production point
  // with the perigee surface. Used both to pick the propagation direction and
  // as the extrapolation point for neutral particles.
  const auto intersection =
      surface
          ->intersect(geoContext, start.position(geoContext), direction,
                      Acts::BoundaryTolerance::Infinite())
          .closest();

  if (start.charge() == 0) {
    // Neutral particle: no helix, linearly extrapolate to the perigee.
    auto lp =
        surface->globalToLocal(geoContext, intersection.position(), direction);
    if (!lp.ok()) {
      return std::nullopt;
    }
    // Direction and momentum are unchanged by a straight-line extrapolation,
    // so only the local position differs from the truth parameters.
    Acts::BoundVector params;
    params << lp.value(), start.phi(), start.theta(), start.qOverP(),
        start.time();
    return Acts::BoundTrackParameters(std::move(surface), params, std::nullopt,
                                      start.particleHypothesis());
  }

  // Charged particle: propagate the truth helix to the perigee surface.
  // BasePropagator::propagateToSurface already forces the perigee surface to be
  // reached; we only steer which side to propagate towards.
  Acts::PropagatorPlainOptions options(geoContext, magFieldContext);
  options.direction =
      Acts::Direction::fromScalarZeroAsPositive(intersection.pathLength());

  auto result = propagator.propagateToSurface(start, *surface, options);
  if (!result.ok()) {
    return std::nullopt;
  }
  return std::move(*result);
}

}  // namespace ActsExamples
