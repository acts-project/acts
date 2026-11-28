// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <functional>
#include <optional>

namespace Acts {

/// Per-call contexts, seeding reference and external fit prior.
///
/// The seeding reference and fit constraint are independent. A reference can
/// guide candidate finding without constraining the subsequent fit:
/// @code{.cpp}
/// Acts::VertexingOptions options(gctx, mctx);
/// options.seedReference = beamSpot;
/// auto vertices = finder.find(tracks, options, state);
/// @endcode
/// To also constrain the fitted positions, set @c options.fitConstraint to the
/// desired vertex prior. This can differ from the seeding reference and can be
/// changed or cleared for each call. A position prior pulls the fitted vertex
/// towards that position; it does not impose a flight-direction constraint.
///
/// @par Migration
/// The former @c constraint member served both purposes. Copy its value to
/// @c seedReference, including its covariance. For a fit using that external
/// prior, also copy it to @c fitConstraint if the former @c useConstraintInFit
/// was true. AMVF's seed-derived mode is described below. Leave @c fitConstraint empty for
/// an unconstrained fit; do not discard a reference merely because the fit is
/// unconstrained. Construct these options with the two contexts, then assign
/// the desired inputs explicitly.
///
/// ZScanVertexFinder::Config::useReferenceForImpactParameters independently
/// controls reference-based seeding. Set it to the old @c useConstraintInFit
/// value to preserve that seeder's treatment of a supplied reference.
/// AdaptiveMultiVertexFinder::Config::useSeedConstraint is now an explicit
/// alternative to an external prior: when migrating that mode, enable it,
/// retain @c seedReference and leave @c fitConstraint empty. The unconstrained
/// mode instead keeps @c useSeedConstraint disabled.
///
/// AdaptiveMultiVertexFitter uses only the contexts from these options; its
/// constraints are specified per candidate in VertexFitProblem.
struct VertexingOptions {
  /// Construct with contexts and no external reference or fit prior.
  ///
  /// @param gctx Geometry context
  /// @param mctx Magnetic field context
  VertexingOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx)
      : geoContext(gctx), magFieldContext(mctx) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// Optional reference for seeding, including its uncertainty. Seeders use
  /// their default behavior when absent. This is not an external fit prior.
  std::optional<Vertex> seedReference;
  /// Optional external position and covariance prior for fitted vertices.
  /// Its covariance must satisfy the fitter's requirements. Seed-only finders
  /// reject this input. AMVF's internal numerical prior is configured
  /// independently and is not represented by this member.
  std::optional<Vertex> fitConstraint;
};

}  // namespace Acts
