// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

template <concepts::point2D point2_t, concepts::point3D point3_t>
struct vertexer;

/// Compute vertices in global frame along the boundary of a surface
///
/// @param ctx geometry context
/// @param sf the surface
/// @param n_seg the number of segments used along arcs
///
/// @returns a vector of vetices (3D points)
template <typename detector_t>
DETRAY_HOST constexpr auto get_global_vertices(
    const typename detector_t::geometry_context &ctx,
    geometry::surface<detector_t> sf, const dindex n_seg) {
  using algebra_t = typename detector_t::algebra_type;
  using point2_t = dpoint2D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;

  auto vertices = sf.template visit_mask<vertexer<point2_t, point3_t>>(n_seg);
  const auto &trf = sf.transform(ctx);

  const std::size_t n_vertices{vertices.size()};
  for (std::size_t i = 0u; i < n_vertices; ++i) {
    vertices[i] = trf.point_to_global(vertices[i]);
  }

  return vertices;
}

/// Generate phi values along an arc
///
/// @param start_phi is the start for the arc generation
/// @param end_phi is the end of the arc generation
/// @param n_seg is the number of segments used to generate the arc
///
/// @return a vector of phi values for the arc
template <concepts::scalar scalar_t>
static inline dvector<scalar_t> phi_values(scalar_t start_phi, scalar_t end_phi,
                                           dindex n_seg) {
  dvector<scalar_t> values;
  values.reserve(n_seg + 1u);
  scalar_t step_phi = (end_phi - start_phi) / static_cast<scalar_t>(n_seg);
  for (unsigned int istep = 0u; istep <= n_seg; ++istep) {
    values.push_back(start_phi + static_cast<scalar_t>(istep) * step_phi);
  }
  return values;
}

/// Create a r-phi polygon from principle parameters
///
/// @param rmin minimum r parameter
/// @param rmax maximum r parameter
/// @param phimin minimum phi parameter
/// @param phimax maximum phi parameters
///
/// @return a polygon representation of the bin
template <concepts::scalar scalar_t, concepts::point2D point2_t>
inline std::vector<point2_t> r_phi_polygon(scalar_t rmin, scalar_t rmax,
                                           scalar_t phimin, scalar_t phimax,
                                           unsigned int n_segments = 1u) {
  std::vector<point2_t> r_phi_poly;
  r_phi_poly.reserve(2u * n_segments + 2u);

  scalar_t cos_min_phi = math::cos(phimin);
  scalar_t sin_min_phi = math::sin(phimin);
  scalar_t cos_max_phi = math::cos(phimax);
  scalar_t sin_max_phi = math::sin(phimax);

  // @TODO add phi generators
  r_phi_poly.push_back({rmin * cos_min_phi, rmin * sin_min_phi});
  r_phi_poly.push_back({rmin * cos_max_phi, rmin * sin_max_phi});
  r_phi_poly.push_back({rmax * cos_max_phi, rmax * sin_max_phi});
  r_phi_poly.push_back({rmax * cos_min_phi, rmax * sin_min_phi});

  return r_phi_poly;
}

/// Functor to produce vertices on a mask collection in a mask tuple container.
template <concepts::point2D point2_t, concepts::point3D point3_t>
struct vertexer {
  /// Specialized method to generate vertices per mask group
  ///
  /// @tparam mask_group_t is the type of the mask collection in a mask cont.
  /// @tparam mask_range_t is the type of the according mask range object
  ///
  /// @param masks is the associated (and split out) mask group
  /// @param range is the range list of masks to be processed
  ///
  /// @return a jagged vector of points of the mask vertices (one per mask)
  template <typename mask_group_t, typename mask_range_t>
  dvector<dvector<point3_t>> operator()(const mask_group_t &masks,
                                        const mask_range_t &range,
                                        unsigned int n_segments = 1) {
    dvector<dvector<point3_t>> mask_vertices = {};
    for (auto i : detray::views::iota(range)) {
      mask_vertices.push_back(masks[i].vertices(n_segments));
    }
    return mask_vertices;
  }
};

}  // namespace detray::detail
