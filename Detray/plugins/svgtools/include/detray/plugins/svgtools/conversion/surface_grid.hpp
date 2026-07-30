// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detail/radius_getter.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/concepts.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/logging.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/conversion/grid.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/proto/grid.hpp"

// System include(s)
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace detray::svgtools::conversion {

namespace detail {

/// A functor to access the bins and get the associated surface indices
struct bin_association_getter {
  template <typename group_t, typename index_t, typename volume_t>
  DETRAY_HOST_DEVICE std::vector<std::vector<std::size_t>> operator()(
      [[maybe_unused]] const group_t& group,
      [[maybe_unused]] const index_t index,
      [[maybe_unused]] const volume_t& vol_desc,
      [[maybe_unused]] const std::array<dindex, 2>& search_window) const {
    using accel_t = typename group_t::value_type;

    if constexpr (concepts::surface_grid<accel_t>) {
      using transform3_t = typename accel_t::local_frame_type::transform3_type;
      using algebra_t = typename accel_t::local_frame_type::algebra_type;
      using scalar_t = typename transform3_t::scalar_type;
      using point2_t = typename transform3_t::point2;

      // The actsvg display only works for 2-dimensional grids
      if constexpr (accel_t::dim != 2u) {
        DETRAY_ERROR_HOST("Only 2D grids can be displayed in actvg");
        return {};
      }

      const accel_t grid = group.at(index);
      const std::size_t n_bins{grid.nbins()};

      std::vector<std::vector<std::size_t>> bin_assoc;
      bin_assoc.reserve(n_bins);

      // Create the bin associations
      auto edges0 = grid.template get_axis<0>().bin_edges();
      auto edges1 = grid.template get_axis<1>().bin_edges();

      // In the svg convention the phi axis has to be the second axis to
      // loop over
      constexpr bool is_cyl{
          std::is_same_v<typename accel_t::local_frame_type,
                         detray::cylindrical2D<algebra_t>> ||
          std::is_same_v<typename accel_t::local_frame_type,
                         detray::concentric_cylindrical2D<algebra_t>>};
      if constexpr (is_cyl) {
        edges0.swap(edges1);
      }

      for (std::size_t i = 1u; i < edges0.size(); ++i) {
        scalar_t p0 = 0.5f * (edges0.at(i) + edges0.at(i - 1));

        for (std::size_t j = 1u; j < edges1.size(); ++j) {
          scalar_t p1 = 0.5f * (edges1.at(j) + edges1.at(j - 1));

          // Create the bin center position estimates for detray
          // (swap cylinder coordinates back)
          point2_t bin_center{p0, p1};
          if constexpr (is_cyl) {
            bin_center = {p1, p0};
          }

          // Get all the bin entries and calculate the loc index
          std::vector<std::size_t> entries;
          for (const auto& sf_desc : grid.search(bin_center, search_window)) {
            // actsvg expects the sensitive surfaces to be numbered
            // starting from zero (per volume)
            entries.push_back(vol_desc.to_local_sf_index(sf_desc.index()));
          }

          // Remove duplicates
          std::ranges::sort(entries);
          auto last = std::unique(entries.begin(), entries.end());
          entries.erase(last, entries.end());

          bin_assoc.push_back(std::move(entries));
        }
      }

      return bin_assoc;
    }

    return {};
  }
};

}  // namespace detail

/// @brief Converts a the grid of a detray volume to a actsvg proto grid.
///
/// @param detector the detector
/// @param index the index of the grid's volume
/// @param view the view
/// @param style the style settings
///
/// @returns a proto grid
template <typename detector_t, typename view_t>
auto surface_grid(const detector_t& detector, const dindex index,
                  const view_t& view,
                  const styling::grid_style& style =
                      styling::tableau_colorblind::grid_style) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;
  using geo_object_ids = typename detector_t::geo_obj_ids;

  const auto vol_desc = detector.volume(index);
  const auto link = vol_desc.template accel_link<geo_object_ids::e_sensitive>();

  // Proactively calculate the reference radius for a cylinder grid
  // (will only be used if the volume actually holds a barrel grid)

  // Get the the radii of the volume portal surfaces
  std::vector<scalar_t> radii{};
  const auto vol = tracking_volume{detector, vol_desc};

  // Passive surfaces could be in the brute force finder, but no
  // sensitive surfaces, since the volume has a grid. Their radii are,
  // however, always within the interval of the portal radii
  for (const auto& pt_desc : vol.portals()) {
    auto r = detector.mask_store()
                 .template visit<detray::detail::outer_radius_getter>(
                     pt_desc.mask());
    if (r.has_value()) {
      radii.push_back(*r);
    }
  }

  scalar_t cyl_ref_radius{0.f};
  if (!radii.empty()) {
    scalar_t inner_r = *std::ranges::min_element(radii);
    scalar_t outer_r = *std::ranges::max_element(radii);

    cyl_ref_radius = 0.5f * static_cast<actsvg::scalar>(inner_r + outer_r);
  }

  return svgtools::conversion::grid(detector.accelerator_store(), link, view,
                                    cyl_ref_radius, style);
}

/// @brief Get the surfaces indices that are registered in the bin neighborhoods
///
/// @param detector the detector
/// @param vol the volume that holds the grid and surfaces
/// @param offset transform a global surface index to a local one for the volume
///
/// @returns a vector of surface indices per neighborhood
template <typename detector_t>
std::vector<std::vector<std::size_t>> get_bin_association(
    const detector_t& det, const detray::tracking_volume<detector_t>& vol,
    const std::array<dindex, 2>& search_window = {2u, 2u}) {
  using geo_object_ids = typename detector_t::geo_obj_ids;

  const auto& vol_desc = det.volume(vol.index());

  if (const auto& link =
          vol_desc.template accel_link<geo_object_ids::e_sensitive>();
      !link.is_invalid()) {
    return det.accelerator_store()
        .template visit<detail::bin_association_getter>(link, vol_desc,
                                                        search_window);
  }

  return {};
}

}  // namespace detray::svgtools::conversion
