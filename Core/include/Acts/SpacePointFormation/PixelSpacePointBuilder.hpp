// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {

class GeometryContext;
class Surface;

namespace PixelSpacePointBuilder {

/// @brief Transform a local covariance to z and r
///
/// Only for surfaces whose reference frame depends on the geometry context
/// alone, i.e. plane and disc, so that a caller can look it up once per module
/// with `Surface::referenceFrame(gctx)`.
///
/// @note A common cheaper shortcut is to split a single variance by the module
/// normal, `varZ = v * (1 - nz²)` and `varR = v * nz²`, dropping the second
/// local variance and the correlation.
///
/// @param rotLocalToGlobal The surface reference frame
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return The (z, r) covariance
SquareMatrix2 computeCovarianceZR(const RotationMatrix3& rotLocalToGlobal,
                                  const Vector3& spacePoint,
                                  const SquareMatrix2& localCov);

/// @brief Get z and r covariance from the local position and covariance
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface The surface associated
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return (z, r) components of the global covariance
/// @deprecated Use computeCovarianceZR with a reference frame from
///             Surface::referenceFrame
[[deprecated(
    "Use computeCovarianceZR with a reference frame from "
    "Surface::referenceFrame")]]
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface,
                          const Vector3& spacePoint,
                          const SquareMatrix2& localCov);

}  // namespace PixelSpacePointBuilder

}  // namespace Acts
