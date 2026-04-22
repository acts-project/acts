// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/bin_fillers.hpp"
#include "detray/builders/detail/radius_getter.hpp"
#include "detray/builders/grid_factory.hpp"
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/navigation/accelerators/spatial_grid.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

namespace detray {

/// @brief Build a grid of a certain shape.
///
/// Decorator class to a volume builder that adds a grid as the volumes
/// geometry accelerator structure.
template <typename detector_t, concepts::grid grid_t,
          typename bin_filler_t = fill_by_pos,
          typename grid_factory_t = grid_factory_type<grid_t>>
class grid_builder : public volume_decorator<detector_t> {
  using link_id_t = typename detector_t::volume_type::object_id;

  // Is the given grid type already an acceleration structure?
  using spatial_grid_t =
      std::conditional_t<concepts::surface_accelerator<grid_t>, grid_t,
                         spatial_grid_impl<grid_t>>;
  using spatial_grid_owning_t = typename spatial_grid_t::template type<true>;

 public:
  using detector_type = detector_t;
  using algebra_type = typename detector_t::algebra_type;
  using value_type = typename detector_type::surface_type;
  using scalar_type = dscalar<algebra_type>;

  /// Decorate a volume with a surface accelerator grid
  DETRAY_HOST
  explicit grid_builder(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : volume_decorator<detector_t>(std::move(vol_builder)) {
    DETRAY_VERBOSE_HOST("Add grid builder to volume: " << this->name());

    // The grid builder provides an acceleration structure to the
    // volume, so don't add sensitive surfaces to the brute force method
    if (this->get_builder()) {
      this->has_accel(true);
    }
  }

  /// Should the passive surfaces be added to the grid ?
  void set_add_passives(bool is_add_passive = false) {
    m_add_passives = is_add_passive;

    DETRAY_VERBOSE_HOST("Add passive surfaces to grid: " << std::boolalpha
                                                         << m_add_passives
                                                         << std::noboolalpha);
  }

  /// Set the surface category this grid should contain (type id in the
  /// accelerator link in the volume descriptor)
  void set_type(std::size_t sf_id) { set_type(static_cast<link_id_t>(sf_id)); }

  /// Set the surface category this grid should contain (type id in the
  /// accelerator link in the volume descriptor)
  void set_type(link_id_t sf_id) {
    // Exclude zero, it is reserved for the brute force method
    assert(static_cast<int>(sf_id) > 0);
    // Make sure the id fits in the volume accelerator link
    assert(sf_id < link_id_t::e_size);

    m_id = sf_id;
  }

  /// Delegate init call depending on @param span type
  template <typename grid_shape_t>
  DETRAY_HOST void init_grid(
      const mask<grid_shape_t, algebra_type> &m,
      const darray<std::size_t, grid_t::dim> &n_bins,
      const std::vector<std::pair<typename grid_t::loc_bin_index, dindex>>
          &bin_capacities = {},
      const darray<std::vector<scalar_type>, grid_t::dim> &ax_bin_edges =
          darray<std::vector<scalar_type>, grid_t::dim>()) {
    DETRAY_VERBOSE_HOST("Initialize surface grid...");

    static_assert(
        std::is_same_v<typename grid_shape_t::template local_frame_type<
                           typename detector_t::algebra_type>,
                       typename grid_t::local_frame_type>,
        "Mask has incorrect shape");

    // Build spatial grid from grid utility and (so far empty) mask
    m_grid = spatial_grid_owning_t{m_factory.template new_grid<grid_t>(
                                       m, n_bins, bin_capacities, ax_bin_edges),
                                   typename spatial_grid_t::mask_type{}};

    DETRAY_VERBOSE_HOST("Created empty grid:\n"
                        << DETRAY_TYPENAME(grid_t) << "\n"
                        << m_grid.axes());
  }

  /// Build the empty grid from axis parameters
  DETRAY_HOST void init_grid(
      const std::vector<scalar_type> &spans,
      const std::vector<std::size_t> &n_bins,
      const std::vector<std::pair<typename grid_t::loc_bin_index, dindex>>
          &bin_capacities = {},
      const std::vector<std::vector<scalar_type>> &ax_bin_edges =
          std::vector<std::vector<scalar_type>>()) {
    DETRAY_VERBOSE_HOST("Initialize surface grid...");

    // Build spatial grid from grid utility and (so far empty) mask
    m_grid =
        spatial_grid_owning_t{m_factory.template new_grid<grid_t>(
                                  spans, n_bins, bin_capacities, ax_bin_edges),
                              typename spatial_grid_t::mask_type{}};

    DETRAY_VERBOSE_HOST("Created empty grid:\n"
                        << DETRAY_TYPENAME(grid_t) << "\n"
                        << m_grid.axes());
  }

  /// Fill grid from existing volume using a bin filling strategy
  /// This can also be called without a volume builder
  template <typename volume_type, typename... Args>
  DETRAY_HOST void fill_grid(
      const detector_t &det, const volume_type &vol,
      const typename detector_t::geometry_context ctx = {},
      const bin_filler_t bin_filler = {}, Args &&...args) {
    bin_filler(m_grid, det, vol, ctx, args...);
  }

  /// Fill grid from externally provided surfaces - temporary solution until
  /// the volume builders can be deployed in the toy detector
  template <typename volume_type, typename surface_container_t,
            typename transform_container_t, typename mask_container_t,
            typename... Args>
  DETRAY_HOST void fill_grid(
      const volume_type &vol, const surface_container_t &surfaces,
      const transform_container_t &transforms, const mask_container_t &masks,
      const typename detector_t::geometry_context ctx = {},
      const bin_filler_t bin_filler = {}, Args &&...args) {
    bin_filler(m_grid, vol, surfaces, transforms, masks, ctx, args...);
  }

  /// Add the volume and the grid to the detector @param det
  DETRAY_HOST
  auto build(detector_t &det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    DETRAY_VERBOSE_HOST("Build surface grid...");

    DETRAY_VERBOSE_HOST(
        " -> Defer to other builders to get complete surface descriptors "
        "first:");

    using surface_desc_t = typename detector_t::surface_type;

    // Add the surfaces (portals and/or passives) that are owned by the vol
    typename detector_t::volume_type *vol_ptr =
        volume_decorator<detector_t>::build(det, ctx);

    DETRAY_VERBOSE_HOST("Resume building with updated surface descriptors");

    // Find the surfaces that should be filled into the grid
    const auto vol = tracking_volume{det, vol_ptr->index()};

    // Grid has not been filled previously, fill it automatically
    if (m_grid.size() == 0u) {
      DETRAY_DEBUG_HOST("-> Filling grid automatically...");

      std::vector<surface_desc_t> surfaces{};
      for (auto &sf_desc : vol.surfaces()) {
        if (sf_desc.is_sensitive() ||
            (m_add_passives && sf_desc.is_passive())) {
          surfaces.push_back(sf_desc);
        }
      }

      this->fill_grid(
          tracking_volume{det, volume_decorator<detector_t>::operator()()},
          surfaces, det.transform_store(), det.mask_store(), ctx);
    } else {
      DETRAY_DEBUG_HOST("-> Grid is prefilled...");

      // The grid is prefilled with surface descriptors that contain the
      // correct LOCAL surface indices per bin (e.g. from file IO).
      // Now add the rest of the linking information, which is only
      // available after the volume builder ran
      for (surface_desc_t &sf_desc : m_grid.all()) {
        assert(!detail::is_invalid_value(sf_desc.index()));

        dindex glob_idx{vol_ptr->to_global_sf_index(sf_desc.index())};
        const auto &new_sf_desc = det.surface(glob_idx);

        assert(new_sf_desc.index() == glob_idx);
        assert(!new_sf_desc.identifier().is_invalid());

        sf_desc = new_sf_desc;
      }
    }

    // For cylinder grids, additional mask information is needed
    if constexpr (concepts::cylindrical<typename grid_t::local_frame_type>) {
      assert(vol_ptr->id() == volume_id::e_cylinder);

      constexpr auto inv_vol_link{
          detray::detail::invalid_value<std::uint8_t>()};
      constexpr auto inf{std::numeric_limits<scalar_type>::max()};

      DETRAY_DEBUG_HOST("Find radius for cylinder grid...");

      // Passive surfaces could be in the brute force finder, but no
      // sensitive surfaces, since the volume has a grid. Their radii are,
      // however, always within the interval of the portal radii
      std::vector<scalar_type> radii{};
      for (auto pt_desc : vol.portals()) {
        auto r = det.mask_store().template visit<detail::outer_radius_getter>(
            pt_desc.mask());

        if (r.has_value()) {
          DETRAY_DEBUG_HOST("Found radius " << *r
                                            << " mm of portal: " << pt_desc);
          radii.push_back(*r);
        }
      }

      scalar_type grid_r{0.f};
      if (!radii.empty()) {
        const scalar_type inner_r{*std::ranges::min_element(radii)};
        const scalar_type outer_r{*std::ranges::max_element(radii)};

        grid_r = 0.5f * (inner_r + outer_r);
      }

      typename spatial_grid_t::mask_type cyl_mask{inv_vol_link, grid_r, -inf,
                                                  inf};

      DETRAY_DEBUG_HOST("Setting mask for cylinder grid: " << cyl_mask);

      m_grid.mask(cyl_mask);
    }

    constexpr auto gid{types::id<typename detector_t::accel, spatial_grid_t>};
    const dindex grid_idx{det.accelerator_store().template size<gid>()};

    DETRAY_DEBUG_HOST("Adding grid to volume in detector. Surface type: "
                      << m_id << ", grid: " << gid
                      << ", grid idx: " << grid_idx);

    // Set: contained surface type, grid type, grid instance index
    vol_ptr->set_accel_link(m_id, gid, grid_idx);

    // Add to detector
    det._accelerators.template push_back<gid>(m_grid);

    DETRAY_DEBUG_HOST("Accelerator link: " << vol_ptr->accel_link());
    DETRAY_VERBOSE_HOST("Successfully built "
                        << gid << " for volume: " << this->name());

    DETRAY_DEBUG_HOST("Finished grid: " << m_grid);

    return vol_ptr;
  }

  /// @returns access to the new grid
  DETRAY_HOST
  auto &get() { return m_grid; }

 private:
  link_id_t m_id{link_id_t::e_sensitive};
  grid_factory_t m_factory{};
  // Data owning grid type, so that surface data can be filled into memory
  spatial_grid_owning_t m_grid{};
  bin_filler_t m_bin_filler{};
  bool m_add_passives{false};
};

/// Grid builder from single components
template <typename detector_t,
          template <class, template <std::size_t> class,
                    typename> class grid_factory_t,
          typename grid_shape_t, typename bin_t,
          template <std::size_t> class serializer_t,
          axis::bounds e_bounds = axis::bounds::e_closed,
          template <typename, typename> class... binning_ts>
using grid_builder_type = grid_builder<
    detector_t,
    typename grid_factory_t<bin_t, serializer_t,
                            typename detector_t::algebra_type>::
        template grid_type<axes<grid_shape_t, e_bounds, binning_ts...>>,
    grid_factory_t<bin_t, serializer_t, typename detector_t::algebra_type>>;

}  // namespace detray
