// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <optional>

namespace Acts::Experimental {

/// Fit of the GBTS tracking filter to a seed, at its innermost space point.
struct GbtsSeedFit final {
  /// Whether the seed has a fit of its own.
  bool valid{false};
  /// Position at the innermost space point.
  Vector3 position{Vector3::Zero()};
  /// Azimuthal angle of the direction.
  double phi{0};
  /// dz/dr of the direction.
  double cotTheta{0};
  /// Signed curvature in the transverse plane.
  double curvature{0};
};

/// Estimate free track parameters from the fit of the GBTS tracking filter.
///
/// The momentum uses the field bending the seed in the transverse plane,
/// Bz - cot(theta) * Br, halfway between its innermost and outermost space
/// points. The fitted curvature holds for the seed as a whole, so the field at
/// its innermost space point would not do where the field varies along it.
///
/// @param fit Fit of the tracking filter
/// @param outerPosition Position of the outermost space point of the seed
/// @param magneticField Magnetic field
/// @param fieldCache Magnetic field cache
/// @return Free track parameters, or nothing without a fit or a field
std::optional<FreeVector> freeParametersFromGbtsSeedFit(
    const GbtsSeedFit& fit, const Vector3& outerPosition,
    const MagneticFieldProvider& magneticField,
    MagneticFieldProvider::Cache& fieldCache);

}  // namespace Acts::Experimental
