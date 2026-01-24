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

/// @brief Get z and r covariance from the local position and covariance
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface The surface associated
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return (z, r) components of the global covariance
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface,
                          const Vector3& spacePoint,
                          const SquareMatrix2& localCov);

}  // namespace PixelSpacePointBuilder

}  // namespace Acts
