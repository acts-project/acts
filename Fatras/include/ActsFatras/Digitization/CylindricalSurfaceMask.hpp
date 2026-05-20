// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsFatras {

/// Clip a 2-D segment that lives in the cylinder readout frame (rPhi, z)
/// to the cylinder bounds, which form an axis-aligned rectangle:
///   rPhi ∈ [R * (averagePhi - halfPhi), R * (averagePhi + halfPhi)]
///   z    ∈ [-halfLengthZ, halfLengthZ]
///
/// Implemented with the Liang–Barsky line-clipping algorithm so that
/// segments crossing one or two boundaries are clipped correctly and
/// fully-outside segments are reported as a masking error.
struct CylindricalSurfaceMask {
  /// Shorthand for a 2-d segment;
  using Segment2D = std::array<Acts::Vector2, 2>;

  /// Apply the cylinder mask to the segment.
  ///
  /// @param surface The CylinderSurface that was hit
  /// @param segment The (rPhi, z) segment from CylindricalSurfaceDrift
  ///
  /// @return The clipped segment or DigitizationError::MaskingError /
  ///         UndefinedSurface on failure.
  Acts::Result<Segment2D> apply(const Acts::Surface& surface,
                                const Segment2D& segment) const;
};

}  // namespace ActsFatras
