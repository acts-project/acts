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

#include <vector>

namespace Acts {

class Surface;

/// @brief An interface for reference point generators
///
/// This is used to generate reference points on surfaces e.g.
/// for filling into grids.
struct IReferenceGenerator {
  virtual ~IReferenceGenerator() = default;

  /// Helper to access reference positions for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference points are to be accessed
  ///
  /// @return a vector of reference points for filling
  virtual const std::vector<Vector3> references(
      const GeometryContext& gctx, const Surface& surface) const = 0;
};

}  // namespace Acts
