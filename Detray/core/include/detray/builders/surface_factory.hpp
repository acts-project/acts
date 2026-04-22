// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/unmasked.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_registry.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <concepts>
#include <exception>
#include <memory>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Generates a number of surfaces for a volume (can be portals, passives
/// or sensitives) and fills them into the containers of a volume builder.
///
/// @tparam detector_t the type of detector the volume belongs to.
/// @tparam mask_shape_t the shape of the surface.
template <typename detector_t, typename mask_shape_t>
class surface_factory : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using volume_link_t = typename detector_t::surface_type::navigation_link;
  using scalar_t = dscalar<algebra_t>;

  // Set individual volume link for portals, but only the mother volume index
  // for other surfaces. The container holds one volume link per portal for
  // all surfaces
  using volume_link_collection = std::vector<std::vector<volume_link_t>>;

 public:
  using detector_type = detector_t;

  /// shorthad for a collection of surface data that can be read by a surface
  /// factory
  using surface_data_t = surface_data<detector_t>;
  using sf_data_collection = std::vector<surface_data_t>;

  /// Empty factory.
  surface_factory() = default;

  /// @returns the current number of surfaces that will be built by this
  /// factory
  DETRAY_HOST
  auto size() const -> dindex override {
    check();
    return static_cast<dindex>(m_bounds.size());
  }

  /// @returns the surface types
  DETRAY_HOST
  auto types() const -> const std::vector<surface_id> & { return m_types; }

  /// @returns the mask boundaries currently held by the factory
  DETRAY_HOST
  auto bounds() const
      -> const std::vector<std::vector<std::vector<scalar_t>>> & {
    return m_bounds;
  }

  /// @returns the transforms currently held by the factory
  DETRAY_HOST
  auto transforms() const
      -> const std::vector<typename detector_t::transform3_type> & {
    return m_transforms;
  }

  /// @returns the volume link(s) currently held by the factory
  DETRAY_HOST
  const auto &volume_links() const { return m_volume_links; }

  /// Add all necessary components to the factory for a single surface
  DETRAY_HOST
  void push_back(surface_data_t &&sf_data) override {
    auto [type, vlinks, index, source, bounds, trf] =
        std::move(sf_data).get_data();

    assert(bounds.front().size() == mask_shape_t::boundaries::e_size);

    m_types.push_back(type);
    m_volume_links.push_back(std::move(vlinks));
    m_indices.push_back(index);
    m_sources.push_back(source);
    m_bounds.push_back(std::move(bounds));
    m_transforms.push_back(trf);
  }

  /// Add all necessary components to the factory from bundled surface
  /// data in @param surface_data .
  DETRAY_HOST
  auto push_back(sf_data_collection &&surface_data) -> void override {
    const auto n_surfaces{static_cast<dindex>(size() + surface_data.size())};

    m_volume_links.reserve(n_surfaces);
    m_indices.reserve(n_surfaces);
    m_sources.reserve(n_surfaces);
    m_bounds.reserve(n_surfaces);
    m_transforms.reserve(n_surfaces);

    // Get per-surface data into detector level container layout
    for (auto &sf_data : surface_data) {
      push_back(std::move(sf_data));
    }
  }

  /// Clear old data
  DETRAY_HOST
  auto clear() -> void override {
    m_types.clear();
    m_volume_links.clear();
    m_indices.clear();
    m_sources.clear();
    m_bounds.clear();
    m_transforms.clear();
  }

  /// Generate the surfaces and add them to given data collections.
  ///
  /// @param volume the volume they will be added to in the detector.
  /// @param surfaces the resulting surface objects.
  /// @param transforms the transforms of the surfaces.
  /// @param masks the masks of the surfaces (all of the same shape).
  /// @param ctx the geometry context.
  ///
  /// @returns index range of inserted surfaces in @param surfaces container
  DETRAY_HOST
  auto operator()([[maybe_unused]] typename detector_t::volume_type &volume,
                  [[maybe_unused]]
                  typename detector_t::surface_lookup_container &surfaces,
                  [[maybe_unused]]
                  typename detector_t::transform_container &transforms,
                  [[maybe_unused]] typename detector_t::mask_container &masks,
                  [[maybe_unused]]
                  typename detector_t::geometry_context ctx = {})
      -> dindex_range override {
    DETRAY_VERBOSE_HOST("Add geometric surfaces...");

    // In case the surfaces container is prefilled with other surfaces
    const auto surfaces_offset{static_cast<dindex>(surfaces.size())};

    DETRAY_VERBOSE_HOST("-> Adding " << size() << " surfaces");

    // Nothing to construct
    if (size() == 0u) {
      return {surfaces_offset, surfaces_offset};
    }

    using mask_t = mask<mask_shape_t, algebra_t, volume_link_t>;

    if constexpr (!types::contains<typename detector_t::masks, mask_t>) {
      std::stringstream err_str{};
      err_str << "Could not find mask type '" << mask_shape_t::name
              << "' in detector";

      DETRAY_FATAL_HOST(err_str.str());
      throw std::invalid_argument(err_str.str());
    } else {
      using surface_t = typename detector_t::surface_type;
      using mask_link_t = typename surface_t::mask_link;
      using material_link_t = typename surface_t::material_link;

      constexpr auto mask_id{types::id<typename detector_t::masks, mask_t>};

      // The material will be added in a later step
      constexpr auto no_material{surface_t::material_id::e_none};

      for (const auto [idx, bounds_per_mask] :
           detray::views::enumerate(m_bounds)) {
        // Append the surfaces relative to the current number of
        // surfaces in the stores
        const dindex sf_idx{detail::is_invalid_value(m_indices[idx])
                                ? dindex_invalid
                                : m_indices[idx]};

        // Add transform
        const dindex trf_idx = this->insert_in_container(
            transforms, m_transforms[idx], sf_idx, ctx);

        // Masks are simply appended, since they are distributed onto
        // multiple containers, their ordering is different from the
        // surfaces
        auto v_links_per_mask = m_volume_links[idx];
        assert(v_links_per_mask.size() == bounds_per_mask.size());
        std::size_t n_masks = bounds_per_mask.size();
        if constexpr (std::is_same_v<mask_shape_t,
                                     unmasked<mask_shape_t::dim>>) {
          masks.template emplace_back<mask_id>(
              empty_context{}, v_links_per_mask.front(),
              detail::invalid_value<scalar_t>());
        } else {
          for (std::size_t i = 0u; i < n_masks; ++i) {
            masks.template emplace_back<mask_id>(
                empty_context{}, bounds_per_mask[i], v_links_per_mask[i]);
          }
        }

        // Add surface with all links set (relative to the given
        // containers)
        mask_link_t mask_link{};
        mask_link.set_id(mask_id);
        const auto mask_idx{
            static_cast<dindex>(masks.template size<mask_id>() - n_masks)};
        if constexpr (concepts::index<typename mask_link_t::index_type>) {
          mask_link.set_index(mask_idx);
        } else {
          mask_link.set_index({mask_idx, static_cast<dindex>(n_masks)});
        }
        // If material is present, it is added in a later step
        material_link_t material_link{no_material, dindex_invalid};

        // Add the surface descriptor at the position given by 'sf_idx'
        this->insert_in_container(surfaces,
                                  {surface_t{trf_idx, mask_link, material_link,
                                             volume.index(), m_types[idx]},
                                   m_sources[idx]},
                                  sf_idx);
      }
    }

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// Run internal consistency check of surface data in the builder
  DETRAY_HOST
  void check() const {
    // This should not happen (need same number and ordering of data)
    assert(m_bounds.size() == m_types.size());
    assert(m_bounds.size() == m_volume_links.size());
    assert(m_bounds.size() == m_indices.size());
    assert(m_bounds.size() == m_sources.size());
    assert(m_bounds.size() == m_transforms.size());
  }

  /// Types of the surface (portal|sensitive|passive)
  std::vector<surface_id> m_types{};
  /// Mask volume link (used for navigation)
  volume_link_collection m_volume_links{};
  /// Indices of surfaces for the placement in container
  /// (counted as "volume-local": 0 to n_sf_in_volume)
  std::vector<dindex> m_indices{};
  /// Source links of surfaces
  std::vector<std::uint64_t> m_sources{};
  /// Mask boundaries of surfaces
  std::vector<std::vector<std::vector<scalar_t>>> m_bounds{};
  /// Transforms of surfaces
  std::vector<typename detector_t::transform3_type> m_transforms{};
};

}  // namespace detray
