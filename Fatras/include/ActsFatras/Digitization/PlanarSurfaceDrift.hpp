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
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <tuple>

namespace ActsFatras {

/// The PlanarSurfaceDrift takes an intersection in the nominal surface and
/// projects the ends into the readout surface, which can be at : -1, 0, 1
///
/// A Lorentz drift angle can be applied.
///
struct PlanarSurfaceDrift {
  /// Shorthand for a 2D segment - drifted segment in 2D
  using Segment2D = std::array<Acts::Vector2, 2>;
  /// Shorthand for a 3D segment  - undrifted segment in 3D
  using Segment3D = std::array<Acts::Vector3, 2>;

  /// Drift the full 3D segment onto a surface 2D readout plane.
  ///
  ///
  /// @param gctx The current Geometry context
  /// @param surface The nominal intersection surface
  /// @param thickness The emulated module/depletion thickness
  /// @param pos The position in global coordinates
  /// @param dir The direction in global coordinates
  /// @param driftdir The drift direction in local (surface) coordinates
  ///
  /// @note A drift direction with no perpendicular component will
  /// result in a segment with no lorentz drift or emulate a 3D pixel
  /// sensor.
  ///
  /// @note The readout is alwayws emulated at the central surface,
  /// as the mask will be deployed there, and the measurement is
  /// presented there/
  ///
  /// @return a tuple of the (drifted) Segment on the readout surface
  /// ( @note without masking ), and original 3D segment (in local 3D frame)
  Acts::Result<std::tuple<Segment2D, Segment3D>> toReadout(
      const Acts::GeometryContext& gctx, const Acts::Surface& surface,
      double thickness, const Acts::Vector3& pos, const Acts::Vector3& dir,
      const Acts::Vector3& driftdir = Acts::Vector3(0., 0., 0.)) const;
};

}  // namespace ActsFatras
