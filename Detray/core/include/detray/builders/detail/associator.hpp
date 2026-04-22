// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/utils/concepts.hpp"

// System include(s).
#include <limits>
#include <vector>

namespace detray::detail {

/// Struct that assigns the center of gravity to a rectangular bin
template <detray::concepts::algebra algebra_t>
struct center_of_gravity_rectangle {
  /// Call operator to the struct, allows to chain several chain operators
  /// together
  ///
  /// @param bin_contour The contour description of the bin -> target
  /// @param surface_contour The contour description of the surface -> test
  ///
  /// @note the bin_contour is asummed to be a rectangle
  ///
  /// @return whether this should be associated
  bool operator()(
      const std::vector<dpoint2D<algebra_t>> &bin_contour,
      const std::vector<dpoint2D<algebra_t>> &surface_contour) const {
    using scalar_t = dscalar<algebra_t>;
    using point2_t = dpoint2D<algebra_t>;

    // Check if centre of gravity is inside bin
    point2_t cgs = {0.f, 0.f};
    for (const auto &svtx : surface_contour) {
      cgs = cgs + svtx;
    }
    cgs = 1.f / static_cast<scalar_t>(surface_contour.size()) * cgs;
    scalar_t min_l0 = std::numeric_limits<scalar_t>::max();
    scalar_t max_l0 = -std::numeric_limits<scalar_t>::max();
    scalar_t min_l1 = std::numeric_limits<scalar_t>::max();
    scalar_t max_l1 = -std::numeric_limits<scalar_t>::max();
    for (const auto &b : bin_contour) {
      min_l0 = math::min(b[0], min_l0);
      max_l0 = math::max(b[0], max_l0);
      min_l1 = math::min(b[1], min_l1);
      max_l1 = math::max(b[1], max_l1);
    }

    if (cgs[0] >= min_l0 && cgs[0] < max_l0 && cgs[1] >= min_l1 &&
        cgs[1] < max_l1) {
      return true;
    }

    return false;
  }
};

/// Check if center of mass is inside a generic polygon bin
template <detray::concepts::algebra algebra_t>
struct center_of_gravity_generic {
  /// Call operator to the struct, allows to chain several chain operators
  /// together
  ///
  /// @param bin_contour The contour description of the bin -> target
  /// @param surface_contour The contour description of the surface -> test
  ///
  /// @return whether this should be associated
  bool operator()(
      const std::vector<dpoint2D<algebra_t>> &bin_contour,
      const std::vector<dpoint2D<algebra_t>> &surface_contour) const {
    using scalar_t = dscalar<algebra_t>;
    using point2_t = dpoint2D<algebra_t>;

    // Check if centre of gravity is inside bin
    point2_t cgs = {0.f, 0.f};
    for (const auto &svtx : surface_contour) {
      cgs = cgs + svtx;
    }
    cgs = 1.f / static_cast<scalar_t>(surface_contour.size()) * cgs;

    std::size_t i = 0u;
    std::size_t j = 0u;
    std::size_t num_points = bin_contour.size();

    bool inside = false;
    for (i = 0u, j = num_points - 1u; i < num_points; j = i++) {
      const auto &pi = bin_contour[i];
      const auto &pj = bin_contour[j];
      if ((((pi[1] <= cgs[1]) && (cgs[1] < pj[1])) ||
           ((pj[1] <= cgs[1]) && (cgs[1] < pi[1]))) &&
          (cgs[0] <
           (pj[0] - pi[0]) * (cgs[1] - pi[1]) / (pj[1] - pi[1]) + pi[0])) {
        inside = !inside;
      }
    }
    return inside;
  }
};

/// Check if the edges of the bin and surface contour overlap
template <detray::concepts::algebra algebra_t>
struct edges_intersect_generic {
  /// Call operator to the struct, allows to chain several chain operators
  /// together
  ///
  /// @param bin_contour The contour description of the bin -> target
  /// @param surface_contour The contour description of the surface -> test
  ///
  /// @return whether this should be associated
  bool operator()(
      const std::vector<dpoint2D<algebra_t>> &bin_contour,
      const std::vector<dpoint2D<algebra_t>> &surface_contour) const {
    using scalar_t = dscalar<algebra_t>;
    using point2_t = dpoint2D<algebra_t>;

    auto intersect = [](const point2_t &pi, const point2_t &pj,
                        const point2_t &pk, const point2_t &pl) {
      scalar_t d =
          (pj[0] - pi[0]) * (pl[1] - pk[1]) - (pj[1] - pi[1]) * (pl[0] - pk[0]);

      if (d != 0.f) {
        double r = ((pi[1] - pk[1]) * (pl[0] - pk[0]) -
                    (pi[0] - pk[0]) * (pl[1] - pk[1])) /
                   d;
        double s = ((pi[1] - pk[1]) * (pj[0] - pi[0]) -
                    (pi[0] - pk[0]) * (pj[1] - pi[1])) /
                   d;
        if (r >= 0.f && r <= 1.f && s >= 0.f && s <= 1.f) {
          return true;
        }
      }
      return false;
    };

    // Loop over bin_contour
    for (std::size_t j = 1u; j <= bin_contour.size(); ++j) {
      std::size_t i = j - 1u;
      std::size_t jc = (j == bin_contour.size()) ? 0u : j;
      const auto &pi = bin_contour[i];
      const auto &pj = bin_contour[jc];
      // Loop over surface_contour
      for (std::size_t k = 1u; k <= surface_contour.size(); ++k) {
        std::size_t l = k - 1u;
        std::size_t kc = (k == surface_contour.size()) ? 0u : k;
        const auto &pl = surface_contour[l];
        const auto &pk = surface_contour[kc];
        if (intersect(pi, pj, pk, pl)) {
          return true;
        }
      }
    }
    return false;
  }
};

}  // namespace detray::detail
