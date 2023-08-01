// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @brief Vertex Finder Options
///
template <typename input_track_t>
struct VertexingOptions {
  /// Default constructor is deleted
  VertexingOptions() = delete;

  /// VertexingOptions with context and beam spot constraint
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param bsConstr The beamspot constraint
  /// @param useConstr Boolean indicating whether vertex constraint should be used
  VertexingOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   const Vertex<input_track_t>& bsConstr,
                   const bool useConstr = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        beamSpot(bsConstr),
        useBeamSpotConstraint(useConstr) {
    if (useBeamSpotConstraint && beamSpot.covariance().determinant() == 0.) {
      throw std::invalid_argument(
          "Vertex constraint covariance matrix must be invertible.");
    }
  }

  /// VertexingOptions with context and without beam spot constraint
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  VertexingOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx)
      : geoContext(gctx), magFieldContext(mctx) {
    beamSpot = Vertex<input_track_t>();
    useBeamSpotConstraint = false;
  }

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// The vertex constraint for the fitting
  Vertex<input_track_t> beamSpot;
  /// Boolean indicating whether we use the vertex constraint
  bool useBeamSpotConstraint;
};

}  // namespace Acts
