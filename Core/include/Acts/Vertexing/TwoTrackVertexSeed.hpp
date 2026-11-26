// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// Geometric seed from two local helical trajectories, without a covariance.
struct TwoTrackVertexSeed {
  /// Midpoint of the two trajectory points selected by the estimator.
  Vector3 position = Vector3::Zero();
  /// Three-dimensional separation of those points (not an exact helix DCA).
  double separation = 0.;
};

/// Estimate a two-track vertex in a uniform field parallel to the z axis.
/// Uses transverse circle intersections, or closest transverse circle points
/// when the circles do not intersect. Selects the intersection with the
/// smallest longitudinal separation, breaking ties by the shortest total
/// transverse arc. Each arc is restricted to half a turn from the supplied
/// track position. In zero field, uses the closest points on two straight lines
/// instead. No beamline, primary-vertex or mass constraint is imposed.
/// @param first First track
/// @param second Second track
/// @param gctx Geometry context
/// @param fieldZ Signed axial field in ACTS units
/// @return Seed, or InvalidInput for invalid/unsupported tracks, or SeedingError
/// for parallel straight tracks or concentric circles with no unique seed.
/// @note Curved tracks must both be charged and have nonzero transverse momentum.
/// This is an initial estimate, not a fit or an exact three-dimensional DCA.
Result<TwoTrackVertexSeed> estimateTwoTrackVertexSeed(
    const BoundTrackParameters& first, const BoundTrackParameters& second,
    const GeometryContext& gctx, double fieldZ);

}  // namespace Acts
