// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// @brief Vertex Fitter Options
///
template <typename input_track_t>
struct VertexFitterOptions {
  /// Default contstructor is deleted
  VertexFitterOptions() = delete;

  /// VertexFitterOptions with context
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param vconstr The pointing contraint to a vertex
  VertexFitterOptions(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> mctx,
      const Vertex<input_track_t>& vconstr = Vertex<input_track_t>())

      : geoContext(gctx), magFieldContext(mctx), vertexConstraint(vconstr) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// The vertex constraint for the fitting
  Vertex<input_track_t> vertexConstraint;
};

}  // namespace Acts
