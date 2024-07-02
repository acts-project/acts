// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <utility>
#include <vector>

namespace Acts {

struct SurfaceBinningMatcher {
  /// The binning tolerance parameters
  using Range = std::pair<double, double>;
  std::vector<Range> tolerances{static_cast<int>(binValues), {0., 0.}};

  SurfaceBinningMatcher() = default;

  SurfaceBinningMatcher(const std::vector<Range>& tolpars)
      : tolerances(tolpars) {}

  /// Check function for surface equivalent
  ///
  /// @param gctx [in] gctx the geometry context for this check
  /// @param bValue the binning value for the binning
  /// @param one first surface for checking
  /// @param other second surface for checking
  bool operator()(const Acts::GeometryContext& gctx, Acts::BinningValue bValue,
                  const Acts::Surface* one, const Acts::Surface* other) const {
    // Fast exit
    if (one == other) {
      return true;
    }

    auto oneExt = one->polyhedronRepresentation(gctx, 1).extent();
    auto otherExt = other->polyhedronRepresentation(gctx, 1).extent();

    double oneMin = oneExt.min(bValue);
    double oneMax = oneExt.max(bValue);

    double otherMin = otherExt.min(bValue);
    double otherMax = otherExt.max(bValue);

    return (std::abs(oneMin - otherMin) <= tolerances[bValue].first &&
            std::abs(oneMax - otherMax) <= tolerances[bValue].second);
  }
};

}  // namespace Acts
