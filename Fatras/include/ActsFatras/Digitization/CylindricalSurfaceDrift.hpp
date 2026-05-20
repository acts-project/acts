// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <array>

namespace ActsFatras {

/// Project a 3D hit segment onto the readout side of a CylinderSurface.
///
/// The cylinder readout uses the 2-D local frame (rPhi, z) where
///   - rPhi = R * phi is the tangential coordinate at the cylinder radius,
///   - z is the coordinate along the cylinder axis.
///
/// Endpoints are returned as Vector2(rPhi, z). The optional Lorentz drift is
/// applied analogously to the planar case: the supplied drift direction is
/// interpreted in the readout-local 3-D frame with axes
/// (tangential, axial, radial), so a non-zero radial component (driftDir.z())
/// selects whether the drift is applied to the entry side (radial > 0) or the
/// exit side (radial < 0). A (0, 0, 0) drift direction means "drift to the
/// central readout plane" (i.e. no in-plane shift, like in the planar case).
struct CylindricalSurfaceDrift {
  /// Shorthand for a 2D segment
  using Segment2D = std::array<Acts::Vector2, 2>;

  /// Drift the full 3D segment onto a cylindrical 2D readout surface.
  ///
  /// @param gctx Geometry context
  /// @param surface The CylinderSurface that was hit
  /// @param thickness Emulated module / depletion thickness (radial)
  /// @param pos Hit position in global coordinates
  /// @param dir Hit direction in global coordinates
  /// @param driftDir Drift direction in readout-local 3-D coordinates
  ///                 (tangential, axial, radial); (0, 0, 0) = no in-plane drift
  ///
  /// @return entry/exit point pair in (rPhi, z)
  Segment2D toReadout(const Acts::GeometryContext& gctx,
                      const Acts::Surface& surface, double thickness,
                      const Acts::Vector3& pos, const Acts::Vector3& dir,
                      const Acts::Vector3& driftDir) const;
};

}  // namespace ActsFatras
