// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/IndexedGridFiller.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>

std::vector<std::size_t> Acts::Experimental::detail::binSequence(
    std::array<std::size_t, 2u> minMaxBins, std::size_t expand,
    std::size_t nBins, Acts::detail::AxisBoundaryType type) {
  // Return vector for iterations
  std::vector<std::size_t> rBins;
  /// Helper method to fill a range
  ///
  /// @param lmin the minimum bin
  /// @param lmax the maximum bin
  auto fill_linear = [&](std::size_t lmin, std::size_t lmax) -> void {
    for (std::size_t b = lmin; b <= lmax; ++b) {
      rBins.push_back(b);
    }
  };
  std::size_t bmin = minMaxBins[0u];
  std::size_t bmax = minMaxBins[1u];

  // Open/Bound cases
  if (type != Acts::detail::AxisBoundaryType::Closed) {
    rBins.reserve(bmax - bmin + 1u + 2 * expand);
    // handle bmin:/max expand it down (for bound, don't fill underflow)
    if (type == Acts::detail::AxisBoundaryType::Bound) {
      bmin = (int(bmin) - int(expand) > 0) ? bmin - expand : 1u;
      bmax = (bmax + expand <= nBins) ? bmax + expand : nBins;
    } else if (type == Acts::detail::AxisBoundaryType::Open) {
      bmin = (int(bmin) - int(expand) >= 0) ? bmin - expand : 0u;
      bmax = (bmax + expand <= nBins + 1u) ? bmax + expand : nBins + 1u;
    }
    fill_linear(bmin, bmax);
  } else {
    // Close case
    std::size_t span = bmax - bmin + 1u + 2 * expand;
    // Safe with respect to the closure point, treat as bound
    if (2 * span < nBins && (bmax + expand <= nBins) &&
        (int(bmin) - int(expand) > 0)) {
      return binSequence({bmin, bmax}, expand, nBins,
                         Acts::detail::AxisBoundaryType::Bound);
    } else if (2 * span < nBins) {
      bmin = int(bmin) - int(expand) > 0 ? bmin - expand : 1u;
      bmax = bmax + expand <= nBins ? bmax + expand : nBins;
      fill_linear(bmin, bmax);
      // deal with expansions over the phi boundary
      if (bmax + expand > nBins) {
        std::size_t overstep = (bmax + expand - nBins);
        fill_linear(1u, overstep);
      }
      if (int(bmin) - int(expand) < 1) {
        std::size_t understep = abs(int(bmin) - int(expand));
        fill_linear(nBins - understep, nBins);
      }
      std::sort(rBins.begin(), rBins.end());
    } else {
      // Jump over the phi boundary
      fill_linear(bmax - expand, nBins);
      fill_linear(1, bmin + expand);
      std::sort(rBins.begin(), rBins.end());
    }
  }
  return rBins;
}
