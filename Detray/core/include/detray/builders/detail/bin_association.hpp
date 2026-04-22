// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detail/associator.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/geometry/detail/vertexer.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/utils/grid/populators.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <vector>

namespace detray::detail {

/// Run the bin association of surfaces (via their contour) to a given 2D grid.
///
/// @param context is the context to win which the association is done
/// @param surfaces a range of detector surfaces
/// @param transforms the transforms that belong to the surfaces
/// @param surface_masks the masks that belong to the surfaces
/// @param grid either a cylinder or disc grid to be filled
/// @param tolerance is the bin_tolerance in the two local coordinates
/// @param absolute_tolerance is an indicator if the tolerance is to be
///        taken absolute or relative
template <typename context_t, typename surface_container_t,
          typename transform_container_t, typename mask_container_t,
          concepts::surface_grid grid_t, concepts::scalar scalar_t>
static inline void bin_association(const context_t & /*context*/,
                                   const surface_container_t &surfaces,
                                   const transform_container_t &transforms,
                                   const mask_container_t &surface_masks,
                                   grid_t &grid,
                                   const darray<scalar_t, 2> &bin_tolerance,
                                   bool absolute_tolerance = true) {
  using algebra_t = typename grid_t::local_frame_type::algebra_type;
  using point2_t = dpoint2D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;

  const auto &axis_0 = grid.template get_axis<0>();
  const auto &axis_1 = grid.template get_axis<1>();

  // Disk type bin association
  if constexpr (std::is_same_v<typename grid_t::local_frame_type,
                               polar2D<algebra_t>>) {
    // Run with two different associators: center of gravity and edge
    // intersection
    center_of_gravity_generic<algebra_t> cgs_assoc;
    edges_intersect_generic<algebra_t> edges_assoc;

    // Loop over all bins and associate the surfaces
    for (unsigned int bin_0 = 0u; bin_0 < axis_0.nbins(); ++bin_0) {
      for (unsigned int bin_1 = 0u; bin_1 < axis_1.nbins(); ++bin_1) {
        auto r_borders = axis_0.bin_edges(bin_0);
        auto phi_borders = axis_1.bin_edges(bin_1);

        scalar_t r_add = absolute_tolerance
                             ? bin_tolerance[0]
                             : bin_tolerance[0] * (r_borders[1] - r_borders[0]);
        scalar_t phi_add =
            absolute_tolerance
                ? bin_tolerance[1]
                : bin_tolerance[1] * (phi_borders[1] - phi_borders[0]);

        // Create a contour for the bin
        std::vector<point2_t> bin_contour =
            detail::r_phi_polygon<scalar_t, point2_t>(
                r_borders[0] - r_add, r_borders[1] + r_add,
                phi_borders[0] - phi_add, phi_borders[1] + phi_add);

        // Run through the surfaces and associate them by contour
        for (auto sf : surfaces) {
          // Add only sensitive surfaces to the grid
          if (sf.is_portal()) {
            continue;
          }

          // Unroll the mask container and generate vertices
          const auto &transform = transforms[sf.transform()];

          auto vertices_per_masks =
              surface_masks
                  .template visit<detail::vertexer<point2_t, point3_t>>(
                      sf.mask());

          // Usually one mask per surface, but design allows - a
          // single association  is sufficient though
          for (auto &vertices : vertices_per_masks) {
            if (!vertices.empty()) {
              // Create a surface contour
              std::vector<point2_t> surface_contour;
              surface_contour.reserve(vertices.size());
              for (const auto &v : vertices) {
                auto vg = transform.point_to_global(v);
                surface_contour.push_back({vg[0], vg[1]});
              }
              // The association has worked
              if (cgs_assoc(bin_contour, surface_contour) ||
                  edges_assoc(bin_contour, surface_contour)) {
                grid.template populate<attach<>>({bin_0, bin_1}, sf);
                break;
              }
            }
          }
        }
      }
    }
  } else if constexpr (std::is_same_v<typename grid_t::local_frame_type,
                                      cylindrical2D<algebra_t>> ||
                       std::is_same_v<typename grid_t::local_frame_type,
                                      concentric_cylindrical2D<algebra_t>>) {
    center_of_gravity_rectangle<algebra_t> cgs_assoc;
    edges_intersect_generic<algebra_t> edges_assoc;

    // Loop over all bins and associate the surfaces
    for (unsigned int bin_0 = 0u; bin_0 < axis_0.nbins(); ++bin_0) {
      for (unsigned int bin_1 = 0u; bin_1 < axis_1.nbins(); ++bin_1) {
        auto z_borders = axis_0.bin_edges(bin_0);
        auto phi_borders = axis_1.bin_edges(bin_1);

        scalar_t z_add = absolute_tolerance
                             ? bin_tolerance[0]
                             : bin_tolerance[0] * (z_borders[1] - z_borders[0]);
        scalar_t phi_add =
            absolute_tolerance
                ? bin_tolerance[1]
                : bin_tolerance[1] * (phi_borders[1] - phi_borders[0]);

        scalar_t z_min = z_borders[0];
        scalar_t z_max = z_borders[1];
        scalar_t phi_min_rep = phi_borders[0];
        scalar_t phi_max_rep = phi_borders[1];

        point2_t p0_bin = {z_min - z_add, phi_min_rep - phi_add};
        point2_t p1_bin = {z_min - z_add, phi_max_rep + phi_add};
        point2_t p2_bin = {z_max + z_add, phi_max_rep + phi_add};
        point2_t p3_bin = {z_max + z_add, phi_min_rep - phi_add};

        std::vector<point2_t> bin_contour = {p0_bin, p1_bin, p2_bin, p3_bin};

        // Loop over the surfaces within a volume
        for (auto sf : surfaces) {
          // Add only sensitive surfaces to the grid
          if (sf.is_portal()) {
            continue;
          }

          // Unroll the mask container and generate vertices
          const auto &transform = transforms.at(sf.transform());

          auto vertices_per_masks =
              surface_masks
                  .template visit<detail::vertexer<point2_t, point3_t>>(
                      sf.mask());

          for (auto &vertices : vertices_per_masks) {
            if (!vertices.empty()) {
              // Create a surface contour
              std::vector<point2_t> surface_contour{};
              surface_contour.reserve(vertices.size());
              scalar_t phi_min = std::numeric_limits<scalar_t>::max();
              scalar_t phi_max = -std::numeric_limits<scalar_t>::max();
              // We potentially need the split vertices
              std::vector<point2_t> s_c_neg{};
              std::vector<point2_t> s_c_pos{};
              scalar_t z_min_neg = std::numeric_limits<scalar_t>::max();
              scalar_t z_max_neg = -std::numeric_limits<scalar_t>::max();
              scalar_t z_min_pos = std::numeric_limits<scalar_t>::max();
              scalar_t z_max_pos = -std::numeric_limits<scalar_t>::max();

              for (const auto &v : vertices) {
                const point3_t vg = transform.point_to_global(v);
                scalar_t phi = math::atan2(vg[1], vg[0]);
                phi_min = math::min(phi, phi_min);
                phi_max = math::max(phi, phi_max);
                surface_contour.push_back({vg[2], phi});
                if (phi < 0.f) {
                  s_c_neg.push_back({vg[2], phi});
                  z_min_neg = math::min(vg[2], z_min_neg);
                  z_max_neg = math::max(vg[2], z_max_neg);
                } else {
                  s_c_pos.push_back({vg[2], phi});
                  z_min_pos = math::min(vg[2], z_min_pos);
                  z_max_pos = math::max(vg[2], z_max_pos);
                }
              }
              // Check for phi wrapping
              std::vector<std::vector<point2_t>> surface_contours{};
              if (phi_max - phi_min > constant<scalar_t>::pi &&
                  phi_max * phi_min < 0.f) {
                s_c_neg.push_back({z_max_neg, -constant<scalar_t>::pi});
                s_c_neg.push_back({z_min_neg, -constant<scalar_t>::pi});
                s_c_pos.push_back({z_max_pos, constant<scalar_t>::pi});
                s_c_pos.push_back({z_min_pos, constant<scalar_t>::pi});
                surface_contours.insert(surface_contours.end(),
                                        {s_c_neg, s_c_pos});
              } else {
                surface_contours.push_back(surface_contour);
              }

              // Check the association (with potential splits)
              bool associated = false;
              for (const auto &s_c : surface_contours) {
                if (cgs_assoc(bin_contour, s_c) ||
                    edges_assoc(bin_contour, s_c)) {
                  associated = true;
                  break;
                }
              }

              // Register if associated
              if (associated) {
                typename grid_t::loc_bin_index mbin{bin_0, bin_1};
                grid.template populate<attach<>>(mbin, sf);
                break;
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace detray::detail
