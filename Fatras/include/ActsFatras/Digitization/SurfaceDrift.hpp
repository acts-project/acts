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
/// A single implementation handles all supported surface types; the readout
/// frame is selected internally from `surface.type()`:
///   - Plane / Disc : the Cartesian local frame (x, y), surface normal = local
///     z. (For discs the polar conversion is done downstream in SurfaceMask /
///     Segmentizer, consistent with the historical behaviour.)
///   - Cylinder     : the unrolled readout frame (rPhi, z), surface normal =
///     radial direction. rPhi = R * phi is the tangential arc length at the
///     cylinder radius.
///
/// In every case the in-plane "x"/"y" coordinates carry the same physical
/// (length) units, so the downstream masking and segmentation are identical.
struct SurfaceDrift {
  /// Shorthand for a 2D segment - drifted segment in 2D readout coordinates
  using Segment2D = std::array<Acts::Vector2, 2>;
  /// Shorthand for a 3D segment - undrifted segment in the local 3D frame
  using Segment3D = std::array<Acts::Vector3, 2>;

  /// Drift the full 3D segment onto the surface 2D readout frame.
  ///
  ///
  /// @param gctx The current Geometry context
  /// @param surface The nominal intersection surface
  /// @param thickness The emulated module/depletion thickness
  /// @param pos The position in global coordinates
  /// @param dir The direction in global coordinates
  /// @param driftDir The drift direction in the readout-local frame
  ///                 (plane/disc: local x, y, normal; cylinder: tangential,
  ///                 axial, radial). A direction with no perpendicular
  ///                 component emulates a 3D pixel sensor / no Lorentz drift.
  ///
  /// @note The readout is always emulated at the central surface,
  /// as the mask will be deployed there, and the measurement is
  /// presented there.
  ///
  /// @return a tuple of the (drifted) Segment2D on the readout surface
  /// ( @note without masking ) and the original undrifted 3D segment, or a
  /// DigitizationError if the track is parallel to the surface.
  Acts::Result<std::tuple<Segment2D, Segment3D>> toReadout(
      const Acts::GeometryContext& gctx, const Acts::Surface& surface,
      double thickness, const Acts::Vector3& pos, const Acts::Vector3& dir,
      const Acts::Vector3& driftDir = Acts::Vector3(0., 0., 0.)) const;
};

}  // namespace ActsFatras
