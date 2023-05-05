// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {

/// A struct to access the center position
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
struct CenterReferenceGenerator {
  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    return {surface.center(gctx)};
  }
};

/// A struct to access reference postions based on bin values
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
struct BinningValueReferenceGenerator {
  /// The binning value
  BinningValue bValue = BinningValue::binValues;

  /// Helper to access a reference postion based on binning value
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    return {surface.binningPosition(gctx, bValue)};
  }
};

/// A struct to access generated vertices from surface polyhedrons
/// These vertices are then used to find the bin boundary box for the
/// indexed grid.
///
/// The grid filling then completes the empty bins in between and
/// expands if necessary.
struct PolyhedronReferenceGenerator {
  /// Also use the barycenter
  bool addBarycenter = true;

  /// The number of segments to approximate (1 means extrema points only)
  unsigned int nSegments = 1;

  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    // Create the return  vector
    std::vector<Vector3> rPositions;
    auto pHedron = surface.polyhedronRepresentation(gctx, nSegments);
    rPositions.insert(rPositions.end(), pHedron.vertices.begin(),
                      pHedron.vertices.end());
    // Add the barycenter if configured
    if (addBarycenter) {
      Vector3 bc(0., 0., 0.);
      std::for_each(rPositions.begin(), rPositions.end(),
                    [&](const auto& p) { bc += p; });
      bc *= 1. / rPositions.size();
      rPositions.push_back(bc);
    }
    return rPositions;
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
